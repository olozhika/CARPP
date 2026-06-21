import numpy as np
from astropy.io import fits
from scipy.ndimage import shift, rotate, zoom
from scipy.optimize import minimize
from .config import default_config
from .models import convolve_gaussian

def read_fits_image(filepath):
    """读取 FITS 文件并返回数据和头文件"""
    with fits.open(filepath) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    return data, header

def write_fits_image(filepath, data, header):
    """将数据写入 FITS 文件"""
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(filepath, overwrite=True)

def get_image_peak(data):
    """获取图像峰值及其坐标"""
    idx = np.unravel_index(np.argmax(data, axis=None), data.shape)
    return data[idx], idx[1], idx[0]  # val, x, y

def align_objective(params, img, ref_img):
    """图像对齐的优化目标函数 (替代 sumsq.pro)"""
    dx, dy = params
    shifted_img = shift(img, [dy, dx], mode='constant', cval=0.0)
    # 计算互相关或最小二乘
    # IDL 中使用的是 -total(as * r)
    return -np.sum(shifted_img * ref_img)

def align_images(img, ref_img):
    """
    通过最小二乘法对齐两个图像。
    对应 IDL 中的 coalign.pro。
    """
    initial_guess = [0.0, 0.0]
    res = minimize(align_objective, initial_guess, args=(img, ref_img), method='Powell')
    dx, dy = res.x
    aligned_img = shift(img, [dy, dx], mode='constant', cval=0.0)
    return aligned_img, [dx, dy]

from astropy.wcs import WCS

def get_beam_from_header(header):
    """
    从 FITS Header 中尝试提取 Beam Size (FWHM, arcsec)。
    优先查找 BMAJ, 然后是 RESOLUTION, RESLN, BEAM 等。
    """
    # 常用关键字 (大小写不敏感处理可通过 header.get 并忽略大小写，但 astropy header 默认处理)
    keys = ['BMAJ', 'FWHM', 'RESOLUTION', 'RESOLU', 'RESLN', 'BEAM']
    
    for key in keys:
        val = header.get(key)
        if val is not None:
            # BMAJ 通常以度为单位
            if key == 'BMAJ':
                return val * 3600.0
            # 其他关键字可能已经是 arcsec 了，这里假设除了 BMAJ 都是 arcsec
            # 如果是度，可能需要进一步判断量级
            if val < 0.1: # 极小值通常对应度
                 return val * 3600.0
            return val
            
    return None

def preprocess_images(image_files, fwhm_beams=None, back_levels=None, n_cells=61, center_type='pixel', center_val=None, cfg=default_config):
    """
    预处理输入图像：根据用户指定的中心坐标进行对齐。
    fwhm_beams: 如果为 None，尝试从 Header 提取。
    """
    processed_images = []
    subpixel_offsets = []
    
    if back_levels is None:
        back_levels = [0.0] * len(image_files)
        
    # 以第一张图为基准确定基准信息
    ref_data, ref_header = read_fits_image(image_files[0])
    ref_wcs = WCS(ref_header)
    
    # 获取像素大小
    if 'CDELT2' in ref_header:
        pixel_size_arcsec = abs(ref_header['CDELT2']) * 3600.0
    elif 'CD2_2' in ref_header:
        pixel_size_arcsec = abs(ref_header['CD2_2']) * 3600.0
    else:
        # 降级方案
        pixel_size_arcsec = 1.0 # 默认 1"
    icent = (n_cells - 1) // 2

    # 统一转换到天球坐标 (RA, Dec) 作为基准中心
    if center_type == 'sky':
        ref_ra, ref_dec = center_val
    else:
        # 如果是 pixel 模式，将第一张图的给定像素坐标转为天球坐标
        # 注意：这里假设用户提供的 center_val 是基于 0 的 [x, y]
        ref_ra, ref_dec = ref_wcs.all_pix2world(center_val[0], center_val[1], 0)
    
    for i, f in enumerate(image_files):
        data, header = read_fits_image(f)
        data = data - back_levels[i]
        wcs = WCS(header)
        
        # 1. 计算该图中的原始像素中心
        px, py = wcs.all_world2pix(ref_ra, ref_dec, 0)

        # 2. 重采样 (如果不同波段像素比例不同)
        curr_pixel_size = abs(header['CDELT2']) * 3600.0
        if abs(curr_pixel_size - pixel_size_arcsec) > 1e-4:
            zoom_factor = curr_pixel_size / pixel_size_arcsec
            data = zoom(data, zoom_factor)
            # 缩放后，像素坐标也需要相应缩放
            px = px * zoom_factor
            py = py * zoom_factor
        
        # 3. 计算离散化的裁剪中心和亚像素偏移
        # 整数坐标用于决定如何裁剪图像
        ix, iy = int(np.round(px)), int(np.round(py))
        
        # dx, dy 是物理中心相对于裁剪图像中心像素 [icent, icent] 的偏差
        dx = px - ix
        dy = py - iy
        
        # 4. 整数坐标裁剪
        y_start, y_end = iy - icent, iy + icent + 1
        x_start, x_end = ix - icent, ix + icent + 1
        
        crop = np.zeros((n_cells, n_cells))
        y_s = max(0, y_start); y_e = min(data.shape[0], y_end)
        x_s = max(0, x_start); x_e = min(data.shape[1], x_end)
        
        crop_part = data[y_s:y_e, x_s:x_e]
        cy_s = max(0, -y_start); cy_e = cy_s + (y_e - y_s)
        cx_s = max(0, -x_start); cx_e = cx_s + (x_e - x_s)
        crop[cy_s:cy_e, cx_s:cx_e] = crop_part
        
        # 5. 单位转换 (Jy/beam -> Jy/pixel)
        header_fwhm = get_beam_from_header(header)
        curr_fwhm = None
        
        if fwhm_beams is not None and fwhm_beams[i] is not None:
            curr_fwhm = fwhm_beams[i]
            if header_fwhm is not None:
                if abs(curr_fwhm - header_fwhm) / header_fwhm > 0.01: # 差异大于1%
                    print(f"Warning (Band {i}): Input FWHM ({curr_fwhm}\") differs from Header FWHM ({header_fwhm:.2f}\"). Using Input value.")
        else:
            curr_fwhm = header_fwhm
            if curr_fwhm is not None:
                print(f"Info (Band {i}): No input FWHM. Automatically detected from header: {curr_fwhm:.2f} arcsec.")
            
        if curr_fwhm is None:
            print(f"Warning (Band {i}): No beam info found for {image_files[i]}, using 1.0 arcsec placeholder.")
            curr_fwhm = 1.0
            
        fwhm_pix = curr_fwhm / pixel_size_arcsec
        beam_area_pix = (np.pi * fwhm_pix**2) / (4.0 * np.log(2.0))
        crop = crop / beam_area_pix
        
        processed_images.append(crop)
        # 记录每张图相对于中心像素的微小偏移
        subpixel_offsets.append((dx, dy))
        
    return np.array(processed_images), subpixel_offsets, pixel_size_arcsec

def read_sed_data(sed_file):
    """
    读取 SED 文本数据。
    对应 IDL 中的 sedread.pro。
    """
    # 这是一个简化的解析器，实际可能需要根据文件格式微调
    wavelengths = []
    fluxes = []
    errors = []
    
    with open(sed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))
                errors.append(float(parts[2]))
                
    return np.array(wavelengths), np.array(fluxes), np.array(errors)
