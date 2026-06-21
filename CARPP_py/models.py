import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
from .config import default_config
from .physics import get_layers_at_b, radiative_transfer

def generate_core_image(temp_prof, dens_prof, rout_pc, beta, lambda_microns, n_layers, nx, pixel_size_pc, offset=(0,0), cfg=default_config):
    """
    生成球形冷核模型的 2D 强度分布图像。
    offset: (dx, dy) 亚像素偏移量
    """
    shells = np.zeros((6, n_layers))
    shells[0, :] = temp_prof
    shells[1, :] = beta
    rout_array = (np.arange(n_layers) + 1.0) * (rout_pc / n_layers)
    shells[2, :] = rout_array
    shells[3, :] = dens_prof
    
    rmax = rout_array[-1]
    lambda_cm = lambda_microns / 1.0e4
    
    # 1D 亮度剖面采样：至少采样到 rmax，步长建议与像素大小相当或更细
    # 为了保证插值质量，我们采样到 rmax
    n_b = int(np.ceil(rmax / pixel_size_pc)) + 2 # 稍微多一点点
    b_array = np.linspace(0, rmax, n_b)
    
    bnu_b = np.zeros(n_b)
    for i in range(n_b):
        layers_data, n_l = get_layers_at_b(shells, n_layers, b_array[i])
        if layers_data is not None:
            bnu_b[i] = radiative_transfer(n_l, layers_data, lambda_cm, cfg)
            
    # --- 亚像素偏移处理 --- (核心修改)
    # dx, dy 是像素坐标偏移。
    dx, dy = offset
    center_pix = (nx - 1) / 2.0
    y, x = np.indices((nx, nx))
    
    # 将模型中心相对于像素网格中心平移 (dx, dy)
    dist_map = np.sqrt((x - (center_pix + dx))**2 + (y - (center_pix + dy))**2) * pixel_size_pc
    
    f_interp = interp1d(b_array, bnu_b, kind='linear', bounds_error=False, fill_value=0.0)
    core_img = f_interp(dist_map)
    
    # 像素立体角应该基于实际像素大小 pixel_size_pc
    omega_pix = (pixel_size_pc * cfg.pc / cfg.distance_cm)**2
    core_img_jy = core_img * omega_pix / cfg.jy
    
    return core_img_jy

def get_core_model(temp_prof, dens_rel_prof, cloud_mass_msun, rout_pc, grain_size_microns, 
                   beta, distance_pc, lambda_microns, n_layers, nx, pixel_size_pc, offset=(0,0), cfg=default_config):
    """
    增加参数: offset=(dx, dy), pixel_size_pc
    """
    cfg.update_derived(grain_size_microns, distance_pc)
    
    # 计算密度
    rlayer_cm = (rout_pc * cfg.pc) / n_layers
    r_cm = np.arange(n_layers) * rlayer_cm
    dr_cm = rlayer_cm
    volumes = (4.0/3.0) * cfg.pi * ((r_cm + dr_cm)**3 - r_cm**3)
    total_mass_rel = np.sum(dens_rel_prof * volumes) * cfg.mh2
    peak_val = (cloud_mass_msun * cfg.msun) / total_mass_rel
    dens_prof = dens_rel_prof * peak_val
    
    img = generate_core_image(temp_prof, dens_prof, rout_pc, beta, lambda_microns, n_layers, nx, pixel_size_pc, offset, cfg)
    
    return img, peak_val

def convolve_gaussian(image, fwhm_pix):
    """
    使用高斯核对图像进行卷积。
    对应 IDL 中的 conga.pro。
    FWHM = 2 * sqrt(2 * ln 2) * sigma approx 2.355 * sigma
    """
    sigma = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return gaussian_filter(image, sigma=sigma, mode='constant', cval=0.0)
