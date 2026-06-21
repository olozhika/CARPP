import numpy as np
import os
import json
import shutil
from astropy.io import fits
from astropy.wcs import WCS
from .config import default_config
from .models import get_core_model, convolve_gaussian

def save_simulation_metadata(output_dir, params, wavelengths, fwhm_beams, pixel_size_arcsec, n_cells, rout_pc, noise_rms, filenames):
    """
    保存模拟的元数据和配置模板
    """
    # 1. 保存结构化的元数据 (JSON)
    metadata = {
        "physical_parameters": {
            "Mass_Msun": float(params[0]),
            "Beta": float(params[1]),
            "Alpha": float(params[2]),
            "T0_K": float(params[3]),
            "T1_K": float(params[4]),
            "rt_pc": float(params[5]),
            "r0_pc": float(params[6]),
            "rout_pc": float(rout_pc)
        },
        "observation_settings": {
            "wavelengths_um": [float(w) for w in wavelengths],
            "fwhm_beams_arcsec": [float(f) for f in fwhm_beams],
            "pixel_size_arcsec": float(pixel_size_arcsec),
            "n_cells": int(n_cells),
            "noise_rms_jy_beam": [float(x) for x in noise_rms] if not np.isscalar(noise_rms) else float(noise_rms),
            "filenames": filenames
        }
    }
    
    with open(os.path.join(output_dir, "simulation_metadata.json"), 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=4)

    # 2. 复制并修改设置模板
    template_path = os.path.join(os.path.dirname(__file__), "carpp_settings_template.txt")
    target_settings = os.path.join(output_dir, "carpp_settings.txt")
    
    if os.path.exists(template_path):
        with open(template_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        new_lines = []
        center_xy = (n_cells - 1) / 2.0
        
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('image_files ='):
                new_lines.append(f"image_files = {filenames} # Updated by Simulator\n")
            elif stripped.startswith('wavelengths_microns ='):
                wave_list = [float(w) for w in wavelengths]
                new_lines.append(f"wavelengths_microns = {wave_list}\n")
            elif stripped.startswith('fwhm_arcsec ='):
                fwhm_list = [float(f) for f in fwhm_beams]
                new_lines.append(f"fwhm_arcsec = {fwhm_list}\n")
            elif stripped.startswith('noise_rms ='):
                noise_val = [float(x) for x in noise_rms] if not np.isscalar(noise_rms) else [float(noise_rms)]*len(wavelengths)
                new_lines.append(f"noise_rms = {noise_val}\n")
            elif stripped.startswith('rout_pc ='):
                new_lines.append(f"rout_pc = {float(rout_pc)}\n")
            elif stripped.startswith('n_cells ='):
                new_lines.append(f"n_cells = {int(n_cells)}\n")
            elif stripped.startswith('center_val ='):
                new_lines.append(f"center_val = [{center_xy}, {center_xy}]\n")
            else:
                new_lines.append(line)
                
        with open(target_settings, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        print(f"Generated fit-ready settings: {target_settings}")

def simulate_observation(params, wavelengths, fwhm_beams, pixel_size_arcsec, n_cells, 
                        rout_pc=None, noise_rms=0.0, output_dir="simulated_data", cfg=default_config):
    """
    生成模拟云核观测图像并保存为 FITS。
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    cloud_mass, beta, alpha, t0, t1, rt, r0 = params
    n_layers = 100 # 生成模型时的解析层数
    
    # 转换噪声为列表
    if np.isscalar(noise_rms):
        noise_rms_list = [noise_rms] * len(wavelengths)
    else:
        noise_rms_list = noise_rms

    distance_pc = cfg.distance_pc
    pixel_size_pc = (pixel_size_arcsec / 206265.0) * distance_pc
    
    # 如果用户没给 rout_pc，默认使用满幅图像半径
    if rout_pc is None:
        rout_pc = ((n_cells - 1) / 2.0) * pixel_size_pc
    
    # 计算密度和温度剖面
    r_layers = (np.arange(n_layers) + 0.5) * (rout_pc / n_layers)
    dens_rel = 1.0 / (1.0 + (r_layers / r0)**alpha)
    temp_prof = t1 + (t0 - t1) / (1.0 + (r_layers / rt)**2.0)

    simulated_files = []
    filenames_only = []
    actual_noise_levels = []

    for i, wave in enumerate(wavelengths):
        fwhm_arcsec = fwhm_beams[i]
        fwhm_pix = fwhm_arcsec / pixel_size_arcsec
        
        # 1. 核心模型生成 (Jy/pixel)
        pad = int(np.ceil(fwhm_pix * 3.0))
        nx_model = n_cells + 2 * pad
        
        model_img_large, _ = get_core_model(
            temp_prof, dens_rel, cloud_mass, rout_pc, 
            cfg.grain_size_cm*1e4, beta, cfg.distance_pc,
            wave, n_layers, nx_model, pixel_size_pc,
            offset=(0,0), cfg=cfg
        )
        
        # 2. 卷积
        model_conv_large = convolve_gaussian(model_img_large, fwhm_pix)
        
        # 3. 裁剪
        c_idx = (nx_model - 1) / 2.0
        icent = (n_cells - 1) / 2.0
        y_s = int(np.round(c_idx - icent))
        x_s = int(np.round(c_idx - icent))
        model_jy_pix = model_conv_large[y_s:y_s+n_cells, x_s:x_s+n_cells]
        
        # 4. 转换回 Jy/beam
        beam_area_pix = (np.pi * fwhm_pix**2) / (4.0 * np.log(2.0))
        model_jy_beam = model_jy_pix * beam_area_pix
        
        # 5. 处理噪声 (支持 "0.05peak" 或 "5%peak")
        current_noise_spec = noise_rms_list[i]
        actual_rms = 0.0
        
        if isinstance(current_noise_spec, str) and "peak" in current_noise_spec.lower():
            peak_val = np.max(model_jy_beam)
            # 提取数字部分
            cleaned = current_noise_spec.lower().replace("peak", "").replace("%", "").strip()
            try:
                factor = float(cleaned)
                if "%" in current_noise_spec:
                    factor /= 100.0
                actual_rms = float(factor * peak_val)
            except ValueError:
                print(f"Warning: Could not parse noise spec '{current_noise_spec}', using 0.")
                actual_rms = 0.0
        else:
            try:
                actual_rms = float(current_noise_spec)
            except (ValueError, TypeError):
                print(f"Warning: Invalid noise level '{current_noise_spec}', using 0.")
                actual_rms = 0.0

        actual_noise_levels.append(actual_rms)
        noise = np.random.normal(0, actual_rms, (n_cells, n_cells))
        obs_sim = model_jy_beam + noise
        
        # 6. 构建 FITS Header
        w = WCS(naxis=2)
        w.wcs.crpix = [n_cells/2.0 + 0.5, n_cells/2.0 + 0.5]
        w.wcs.cdelt = [-pixel_size_arcsec/3600.0, pixel_size_arcsec/3600.0]
        w.wcs.crval = [0.0, 0.0]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        
        header = w.to_header()
        header['BUNIT'] = 'Jy/beam'
        header['BMAJ'] = fwhm_arcsec / 3600.0
        header['BMIN'] = fwhm_arcsec / 3600.0
        header['BPA'] = 0.0
        header['WAVELEN'] = (wave, 'microns')
        header['INSTRUME'] = 'CARPP_Simulator'
        
        fname = f"image-{int(wave)}um.fits"
        filepath = os.path.join(output_dir, fname)
        
        hdu = fits.PrimaryHDU(data=obs_sim.astype(np.float32), header=header)
        hdu.writeto(filepath, overwrite=True)
        simulated_files.append(filepath)
        filenames_only.append(fname)
        
        print(f"Generated simulated image: {filepath} (Noise RMS: {actual_rms:.2e} Jy/beam)")

    # 保存额外信息
    save_simulation_metadata(output_dir, params, wavelengths, fwhm_beams, pixel_size_arcsec, n_cells, rout_pc, actual_noise_levels, filenames_only)

    return simulated_files

def calculate_optical_depths(metadata_or_path, cfg=default_config):
    """
    计算指定仿真配置下的径向层光学厚度和视线（projected）方向的光学厚度分布。
    
    参数:
    metadata_or_path: simulation_metadata.json 的文件路径，或者直接是解析出来的 dict 对象。
    cfg: 配置对象，默认使用全局 default_config
    
    返回:
     一个字典，包含物理尺度、每一层密度和各个波长下的光学厚度。
    """
    if isinstance(metadata_or_path, str):
        with open(metadata_or_path, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
    else:
        metadata = metadata_or_path
        
    p = metadata["physical_parameters"]
    obs = metadata["observation_settings"]
    
    mass = p["Mass_Msun"]
    beta = p["Beta"]
    alpha = p["Alpha"]
    t0 = p["T0_K"]
    t1 = p["T1_K"]
    rt = p["rt_pc"]
    r0 = p["r0_pc"]
    rout_pc = p["rout_pc"]
    
    wavelengths = obs["wavelengths_um"]
    
    n_layers = 100
    r_layers = (np.arange(n_layers) + 0.5) * (rout_pc / n_layers)
    dr_pc = rout_pc / n_layers
    dr_cm = dr_pc * cfg.pc
    
    # 重构气体分子数密度 (H2 cm^-3)
    dens_rel = 1.0 / (1.0 + (r_layers / r0)**alpha)
    
    r_cm = np.arange(n_layers) * dr_cm
    volumes = (4.0/3.0) * cfg.pi * ((r_cm + dr_cm)**3 - r_cm**3)
    total_mass_rel = np.sum(dens_rel * volumes) * cfg.mh2
    peak_val = (mass * cfg.msun) / total_mass_rel
    dens_prof = dens_rel * peak_val  # Physical density in cm^-3
    
    from .physics import get_layers_at_b
    rout_array = (np.arange(n_layers) + 1.0) * dr_pc
    
    # 构造 shells 数据以便调用 get_layers_at_b
    shells = np.zeros((6, n_layers))
    shells[0, :] = t1 + (t0 - t1) / (1.0 + (r_layers / rt)**2.0)
    shells[1, :] = beta
    shells[2, :] = rout_array
    shells[3, :] = dens_prof
    
    results = {
        "wavelengths_um": wavelengths,
        "r_pc": r_layers.tolist(),
        "density_cm3": dens_prof.tolist(),
        "optical_depths": {}
    }
    
    for wave in wavelengths:
        lambda_cm = wave / 1e4
        q_layer = cfg.q350 * (lambda_cm / 350.0e-4)**beta
        # 尘埃不透明度折合成气体的等效不透明度
        opacity_coeff = (cfg.mh2 / cfg.gdr) * (3.0 * q_layer / (4.0 * cfg.grain_size_cm * cfg.rho_d))
        
        # 1. 径向单层增量光学厚度
        radial_inc = dens_prof * opacity_coeff * dr_cm
        
        # 2. 从外边缘到内部累积的径向光学厚度
        radial_cum = np.cumsum(radial_inc[::-1])[::-1]
        
        # 3. 视线方向光学厚度 profile tau_LOS(b)
        los_profile = np.zeros(n_layers)
        for k in range(n_layers):
            b_val = r_layers[k]
            layers_data, n_l = get_layers_at_b(shells, n_layers, b_val)
            if layers_data is not None:
                tau_los_sum = 0.0
                for i in range(n_l):
                    nh2_l = layers_data[3, i]
                    l_pc = layers_data[2, i]
                    l_cm = l_pc * cfg.pc
                    tau_los_sum += nh2_l * opacity_coeff * l_cm
                los_profile[k] = tau_los_sum
                
        # 视线最中心处的总光学厚度 tau_LOS(0)
        layers_data_center, n_l_center = get_layers_at_b(shells, n_layers, 0.0)
        max_tau_los = 0.0
        if layers_data_center is not None:
            for i in range(n_l_center):
                nh2_l = layers_data_center[3, i]
                l_pc = layers_data_center[2, i]
                l_cm = l_pc * cfg.pc
                max_tau_los += nh2_l * opacity_coeff * l_cm
        else:
            max_tau_los = 2.0 * np.sum(radial_inc)
            
        results["optical_depths"][wave] = {
            "radial_incremental": radial_inc.tolist(),
            "radial_cumulative": radial_cum.tolist(),
            "los_profile": los_profile.tolist(),
            "max_tau_los": float(max_tau_los)
        }
        
    return results

