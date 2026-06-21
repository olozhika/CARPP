import os
import time
import numpy as np
from .config import Config
from .preprocessing import preprocess_images
from .fitter import CoreFitter
from .physics import calculate_central_density

def detect_and_parse_sed_file(config_dir, settings=None):
    """
    自动检测并解析配置目录中的辅助积分 SED 文件 / 光谱数据。
    支持：
    1. 传统的以 .txt 结尾且首行包含 '# Spatially integrated SED' 的高积分 SED 文件。
    2. 新的以 .dat / .txt 结尾且包含 '# source = ' 或 '# aperture = ' 等信息的 Spitzer 类似多孔径光谱数据 (e.g. input_spectra.dat)。
    """
    import re
    if not os.path.exists(config_dir):
        return None
    
    for filename in os.listdir(config_dir):
        if filename.endswith(".txt") or filename.endswith(".dat"):
            filepath = os.path.join(config_dir, filename)
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    first_lines = [f.readline() for _ in range(5)]
                
                # Format 1: Old '# Spatially integrated SED' format
                if any("# Spatially integrated SED" in line for line in first_lines if line):
                    width_arcsec = None
                    width_pc = None
                    wavelengths = []
                    fluxes = []
                    reading_data = False
                    
                    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                        for line in f:
                            stripped = line.strip()
                            if not stripped:
                                continue
                            if stripped.startswith('#'):
                                if reading_data and ('seds for' in stripped.lower() or 'observational' in stripped.lower() or '===' in stripped):
                                    break
                                if 'width' in stripped.lower():
                                    match_arcsec = re.search(r'width\s*=\s*([0-9.-eE]+)', stripped, re.IGNORECASE)
                                    if match_arcsec:
                                        width_arcsec = float(match_arcsec.group(1))
                                    match_pc = re.search(r'\(\s*([0-9.-eE]+)\s*pc\)', stripped, re.IGNORECASE)
                                    if match_pc:
                                        width_pc = float(match_pc.group(1))
                                continue
                            
                            parts = stripped.split()
                            if len(parts) >= 2:
                                try:
                                    lam = float(parts[0])
                                    flx = float(parts[1])
                                    wavelengths.append(lam)
                                    fluxes.append(flx)
                                    reading_data = True
                                except ValueError:
                                    if reading_data:
                                        break
                    if len(wavelengths) > 0:
                        return {
                            'filepath': filepath,
                            'width_arcsec': width_arcsec,
                            'width_pc': width_pc,
                            'wavelengths': wavelengths,
                            'fluxes': fluxes,
                            'errors': None
                        }
                
                # Format 2: New multi-aperture dat/txt format (e.g. input_spectra.dat)
                is_spectra_format = False
                for line in first_lines:
                    if line and ("# source" in line.lower() or "# aperture" in line.lower() or "# cdelt" in line.lower()):
                        is_spectra_format = True
                        break
                
                if is_spectra_format:
                    aperture_val = None
                    cdelt_val = None
                    ra_coords = {}
                    dec_coords = {}
                    header_cols = []
                    wavelengths = []
                    reading_data = False
                    
                    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                        for line in f:
                            stripped = line.strip()
                            if not stripped:
                                continue
                            if stripped.startswith('#'):
                                # parse constants like aperture and cdelt
                                if 'aperture' in stripped.lower() and '=' in stripped:
                                    m = re.search(r'aperture\s*=\s*([0-9.-eE]+)', stripped, re.IGNORECASE)
                                    if m:
                                        aperture_val = float(m.group(1))
                                if 'cdelt' in stripped.lower() and '=' in stripped:
                                    m = re.search(r'cdelt\s*=\s*([0-9.-eE]+)', stripped, re.IGNORECASE)
                                    if m:
                                        cdelt_val = float(m.group(1))
                                
                                # parse RA at centers of apertures
                                m_ra = re.search(r'ra_([0-9]+)\s*=\s*([0-9.-eE]+)', stripped, re.IGNORECASE)
                                if m_ra:
                                    idx = int(m_ra.group(1))
                                    ra_coords[idx] = float(m_ra.group(2))
                                
                                # parse DEC at centers of apertures
                                m_dec = re.search(r'dec_([0-9]+)\s*=\s*([0-9.-eE]+)', stripped, re.IGNORECASE)
                                if m_dec:
                                    idx = int(m_dec.group(1))
                                    dec_coords[idx] = float(m_dec.group(2))
                                
                                # parse column headers
                                if 'wave' in stripped.lower() and 'flux_1' in stripped.lower():
                                    header_cols = stripped.lstrip('#').strip().split()
                                continue
                            
                            # Data rows
                            parts = stripped.split()
                            if len(parts) >= 3:
                                # We need to know which aperture is of interest. Since it requires matching the core center,
                                # we process raw parts and we'll resolve the closest aperture indices at the end.
                                try:
                                    row_vals = [float(p) for p in parts]
                                    wavelengths.append(row_vals)
                                    reading_data = True
                                except ValueError:
                                    pass
                    
                    if len(wavelengths) > 0 and reading_data:
                        # Combine RA and Dec to coordinate dict
                        aperture_coords = {}
                        for idx in ra_coords:
                            if idx in dec_coords:
                                aperture_coords[idx] = (ra_coords[idx], dec_coords[idx])
                        
                        # Find core's center RA and Dec from settings
                        ref_ra = None
                        ref_dec = None
                        if settings is not None:
                            center_type = settings.get('center_type', 'pixel')
                            center_val = settings.get('center_val')
                            image_files = settings.get('image_files')
                            if center_type == 'sky':
                                ref_ra, ref_dec = center_val
                            elif center_type == 'pixel' and image_files is not None and len(image_files) > 0:
                                try:
                                    fits_path = image_files[0]
                                    if not os.path.isabs(fits_path):
                                        fits_path = os.path.join(config_dir, fits_path)
                                    if os.path.exists(fits_path):
                                        from astropy.io import fits
                                        from astropy.wcs import WCS
                                        with fits.open(fits_path) as hdul:
                                            header = hdul[0].header
                                            wcs = WCS(header)
                                            ref_ra, ref_dec = wcs.all_pix2world(center_val[0], center_val[1], 0)
                                            ref_ra = float(ref_ra)
                                            ref_dec = float(ref_dec)
                                except Exception as e:
                                    print(f"Warning in matching coordinates for auxiliary spectra: {e}")
                        
                        # Find closest aperture index
                        best_i = None
                        if ref_ra is not None and ref_dec is not None and len(aperture_coords) > 0:
                            min_dist = float('inf')
                            for idx, (ara, adec) in aperture_coords.items():
                                dist = (ara - ref_ra)**2 + (adec - ref_dec)**2
                                if dist < min_dist:
                                    min_dist = dist
                                    best_i = idx
                        
                        if best_i is None:
                            if len(aperture_coords) > 0:
                                sorted_idxs = sorted(list(aperture_coords.keys()))
                                best_i = sorted_idxs[len(sorted_idxs) // 2]
                            else:
                                best_i = 3 # fallback to central aperture
                        
                        print(f"Auto-selected aperture index from spectral file: {best_i}")
                        
                        # Find indices of wave, flux and dflux columns
                        wave_idx = 0
                        flux_idx = 2 * best_i - 1
                        dflux_idx = 2 * best_i
                        
                        if len(header_cols) > 0:
                            try:
                                wave_col_name = 'wave'
                                for c in header_cols:
                                    if c.lower() == 'wave' or c.lower() == 'wavelength':
                                        wave_col_name = c
                                        break
                                wave_idx = header_cols.index(wave_col_name)
                                flux_idx = header_cols.index(f"flux_{best_i}")
                                dflux_idx = header_cols.index(f"dflux_{best_i}")
                            except Exception:
                                pass
                        
                        # Extract final arrays
                        final_lambdas = []
                        final_fluxes = []
                        final_errors = []
                        for row in wavelengths:
                            if len(row) > max(wave_idx, flux_idx, dflux_idx):
                                final_lambdas.append(row[wave_idx])
                                final_fluxes.append(row[flux_idx])
                                final_errors.append(row[dflux_idx])
                        
                        # Calculate aperture width
                        width_arcsec = None
                        if aperture_val is not None and cdelt_val is not None:
                            width_arcsec = aperture_val * cdelt_val
                        
                        width_pc = None
                        if width_arcsec is not None and settings is not None:
                            distance_pc = float(settings.get('distance_pc', 414.0))
                            width_pc = (width_arcsec / 206265.0) * distance_pc
                        
                        if len(final_lambdas) > 0:
                            return {
                                'filepath': filepath,
                                'width_arcsec': width_arcsec,
                                'width_pc': width_pc,
                                'wavelengths': final_lambdas,
                                'fluxes': final_fluxes,
                                'errors': final_errors
                            }
            except Exception as e:
                print(f"Error parsing file {filename} as SED file: {e}")
                pass
    return None

def run_carpp(config_path, method=None, workers=None, bounds=None, initial_guess=None):
    """
    主运行函数：读取配置文件并执行拟合。
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    print(f"--- Loading configuration from: {config_path} ---")
    settings = Config.from_file(config_path)
    
    # 自动搜索并载入伴随高积分 SED 辅助数据
    config_dir = os.path.dirname(os.path.abspath(config_path))
    sed_data = detect_and_parse_sed_file(config_dir, settings=settings)
    if sed_data is not None:
        # 支持限制 SED 波长范围的过滤设置
        sed_min_lambda = float(settings.get('sed_min_lambda', 0.0))
        sed_max_lambda = float(settings.get('sed_max_lambda', 1e6))
        
        keep_idxs = [i for i, lam in enumerate(sed_data['wavelengths']) if sed_min_lambda <= lam <= sed_max_lambda]
        if len(keep_idxs) < len(sed_data['wavelengths']):
            print(f"Filtering supplementary SED wavelengths to range [{sed_min_lambda}, {sed_max_lambda}] microns...")
            sed_data['wavelengths'] = [sed_data['wavelengths'][i] for i in keep_idxs]
            sed_data['fluxes'] = [sed_data['fluxes'][i] for i in keep_idxs]
            if 'errors' in sed_data and sed_data['errors'] is not None:
                sed_data['errors'] = [sed_data['errors'][i] for i in keep_idxs]
        
        print(f"Loaded auxiliary SED file successfully: {os.path.basename(sed_data['filepath'])}")
        print(f"  Wavelengths count: {len(sed_data['wavelengths'])} ({min(sed_data['wavelengths']):.2f} to {max(sed_data['wavelengths']):.2f} um)")
        if sed_data['width_arcsec'] is not None:
            print(f"  Aperture Width   : {sed_data['width_arcsec']:.2f} arcsec")
    
    # 1. 物理环境初始化
    cfg = Config(
        grain_size_microns=settings.get('grain_radius_microns', 0.1),
        distance_pc=settings.get('distance_pc', 400.0)
    )
    
    # 2. 图像预处理 (使用指定的坐标中心)
    print("Preprocessing images (Resampling & Alignment)...")
    n_cells = int(settings['n_cells'])
    default_center = [(n_cells - 1) / 2.0, (n_cells - 1) / 2.0]
    
    obs_images, subpixel_offsets, pixel_size_arcsec = preprocess_images(
        settings['image_files'], 
        settings['fwhm_arcsec'], 
        settings['back_levels'], 
        n_cells,
        center_type=settings.get('center_type', 'pixel'),
        center_val=settings.get('center_val', default_center),
        cfg=cfg
    )
    
    # 3. 拟合引擎准备
    qz1 = settings.get('qz1', 2.0)
    noise_rms = settings.get('noise_rms', [0.1]*len(settings['image_files']))
    cal_error = settings.get('cal_error', 0.1)
    
    # 根据 Header 里的像素比例以及用户距离，计算物理像素大小
    distance_pc = cfg.distance_pc
    pixel_size_pc = (pixel_size_arcsec / 206265.0) * distance_pc
    
    # 智能处理 rout_pc
    rout_pc = settings.get('rout_pc')
    if rout_pc is None or str(rout_pc).lower() == 'none':
        # 如果未指定，则假设云核刚好填满裁剪区域，半径为 (n_cells - 1) / 2
        rout_pc = ((n_cells - 1) / 2.0) * pixel_size_pc
        print(f"rout_pc set to None, automatically calculated as {rout_pc:.4f} pc")
    else:
        rout_pc = float(rout_pc)
    
    print(f"Weighting factor qz1 = {qz1}")
    print(f"Calibration error factor = {cal_error}")
    print(f"Noise RMS levels = {noise_rms}")
    print(f"Cloud outer radius (rout_pc) = {rout_pc:.4f} pc")
    print(f"Calculated pixel size = {pixel_size_pc:.6f} pc (based on FITS scale {pixel_size_arcsec:.3f}\")")
    
    sed_noise_rms = float(settings.get('sed_noise_rms', 0.0))
    sed_cal_error = float(settings.get('sed_cal_error', cal_error)) # 默认使用和图像一样的校准误差
    
    fitter = CoreFitter(obs_images, settings['wavelengths_microns'], settings['fwhm_arcsec'], 
                        subpixel_offsets=subpixel_offsets, qz1=qz1, 
                        noise_rms=noise_rms, cal_error=cal_error, 
                        rout_pc=rout_pc, pixel_size_pc=pixel_size_pc, cfg=cfg,
                        sed_data=sed_data, sed_noise_rms=sed_noise_rms, sed_cal_error=sed_cal_error)
    
    # 4. 执行拟合 (带计时器)
    fit_method = method.lower() if method else settings.get('method', 'powell').lower()
    num_workers = workers if workers is not None else int(settings.get('workers', 1))
    
    # 智能控制参数
    popsize = int(settings.get('popsize', 15))
    maxiter = int(settings.get('maxiter', 100))
    tol = float(settings.get('tol', 0.01))
    
    # MCMC 参数
    nwalkers = int(settings.get('nwalkers', 32))
    nsteps = int(settings.get('nsteps', 1000))
    
    if fit_method != 'powell':
        import multiprocessing
        actual_workers = multiprocessing.cpu_count() if num_workers == -1 else num_workers
        print(f"Parallel/MCMC mode enabled: using {actual_workers} workers.")
    
    fit_bounds = bounds if bounds is not None else settings.get('bounds', None)
    initial_p = initial_guess if initial_guess is not None else settings.get('initial_guess', None)
    if isinstance(initial_p, dict) and 'params' in initial_p:
        initial_p = initial_p['params']
    
    start_time = time.time()  # --- 计时开始 ---
    
    best_fit_params, final_chi2, nit = fitter.fit(
        initial_guess=initial_p, 
        bounds=fit_bounds, 
        method=fit_method, 
        workers=num_workers,
        popsize=popsize,
        maxiter=maxiter,
        tol=tol,
        nwalkers=nwalkers,
        nsteps=nsteps
    )
    
    end_time = time.time()    # --- 计时结束 ---
    elapsed_time = end_time - start_time
    
    # 5. 输出结果
    print("\n" + "="*40)
    print("         CARPP FIT SUMMARY")
    print("="*40)
    param_names = ["Cloud Mass", "Beta", "Alpha", "T0", "T1", "rt", "r0"]
    for name, val in zip(param_names, best_fit_params):
        print(f"{name:15}: {val:.4f}")
    
    print("-" * 40)
    print(f"Reduced Chi-sq : {final_chi2:.4e}")
    print(f"Total Steps    : {nit}")
    print(f"Time Elapsed   : {elapsed_time:.2f} seconds")  # 打印时间
    print(f"Method Used    : {fit_method.upper()}")
    
    # 6. 物理衍生量计算 (中心体密度 nc)
    best_mass = best_fit_params[0]
    best_alpha = best_fit_params[2]
    best_r0 = best_fit_params[6]
    nc_val = calculate_central_density(best_mass, best_alpha, best_r0, rout_pc, cfg=cfg)
    print(f"Central n_H2   : {nc_val:.4e} cm^-3")
    print("="*40)
    
    # 返回一个包含丰富信息的字典
    return {
        'params': best_fit_params,
        'central_density_cm3': nc_val,
        'chi2': final_chi2,
        'steps': nit,
        'time': elapsed_time,
        'method': fit_method,
        'units': {
            'mass': 'M_sun',
            'temp': 'K',
            'scale': 'pc',
            'time': 'sec'
        }
    }

if __name__ == "__main__":
    # 默认运行 data 文件夹下的配置
    run_carpp('data/carpp_settings.txt')
