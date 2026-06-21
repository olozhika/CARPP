import numpy as np
from numba import njit
from .config import default_config

@njit
def planck_function_numba(lambda_cm, temp_k, h, c, k):
    """Numba 加速的普朗克函数"""
    temp_k = max(temp_k, 0.1)
    exponent = (h * c) / (k * lambda_cm * temp_k)
    exponent = min(exponent, 700.0)
    exponent = max(exponent, 1e-10)
    
    term1 = (2.0 * h * c) / (lambda_cm**3.0)
    term2 = 1.0 / (np.exp(exponent) - 1.0 + 1e-20)
    return term1 * term2

@njit
def dust_opacity_numba(lambda_cm, beta, q350):
    """Numba 加速的尘埃不透明度"""
    return q350 * (lambda_cm / 350.0e-4)**beta

@njit
def radiative_transfer_numba(n_layers, layers_data, lambda_cm, pc, h, c, k, q350, mh2, gdr, grain_size, rho_d):
    """
    Numba 加速的辐射传输循环。
    """
    snu_layer = 0.0
    for i in range(n_layers):
        td = layers_data[0, i]
        beta = layers_data[1, i]
        l_layer_cm = layers_data[2, i] * pc
        nh2 = layers_data[3, i]
        
        q_layer = dust_opacity_numba(lambda_cm, beta, q350)
        tau_layer = (nh2 * mh2 / gdr) * (3.0 * q_layer / (4.0 * grain_size * rho_d)) * l_layer_cm
        
        source_b = planck_function_numba(lambda_cm, td, h, c, k)
        
        exp_tau = np.exp(-max(tau_layer, 0.0))
        snu_layer = snu_layer * exp_tau + source_b * (1.0 - exp_tau)
    return snu_layer

def planck_function(lambda_cm, temp_k, cfg=default_config):
    return planck_function_numba(lambda_cm, temp_k, cfg.h, cfg.c, cfg.k)

def dust_opacity(lambda_cm, beta, cfg=default_config):
    return dust_opacity_numba(lambda_cm, beta, cfg.q350)

def radiative_transfer(n_layers, layers_data, lambda_cm, cfg=default_config):
    return radiative_transfer_numba(
        n_layers, layers_data, lambda_cm, 
        cfg.pc, cfg.h, cfg.c, cfg.k, cfg.q350, 
        cfg.mh2, cfg.gdr, cfg.grain_size_cm, cfg.rho_d
    )

def get_layers_at_b(shells, n_shells, b_pc):
    """
    计算给定冲击参数 b 时，视线穿过的所有壳层及其路径长度。
    对应 IDL 中的 layers_b.pro。
    
    shells 结构 [6, n_shells]:
      [0, :] = T
      [1, :] = beta
      [2, :] = 壳层外径 rout [pc]
      [3, :] = 密度 nh2
    """
    if b_pc > shells[2, n_shells-1]:
        # 冲击参数超过了云半径
        return None, 0

    # 寻找视线穿过的最内层索引
    # IDL: while (b gt shells[2,i_center] ) do i_center = i_center+1
    i_center = np.searchsorted(shells[2, :], b_pc)
    
    # 总层数：视线穿过球体，会经过两次（进入和离开），除了最中心
    n_layer = (n_shells - i_center - 1) * 2 + 1
    layers = np.zeros((6, n_layer))
    
    r_center = shells[2, i_center]
    
    if n_layer == 1:
        # 只穿过最内层
        layers[:, 0] = shells[:, i_center]
        layers[2, 0] = 2.0 * np.sqrt(r_center**2 - b_pc**2)
    else:
        # 1. 后半球 (Back half)
        for i in range((n_layer - 1) // 2):
            i_out = n_shells - i - 1
            i_in = i_out - 1
            r_out = shells[2, i_out]
            r_in = shells[2, i_in]
            layers[:, i] = shells[:, i_out]
            # 路径长度计算
            layers[2, i] = np.sqrt(r_out**2 - b_pc**2) - np.sqrt(r_in**2 - b_pc**2)
            
        # 2. 中间层 (Middle shell)
        mid_idx = (n_layer - 1) // 2
        layers[:, mid_idx] = shells[:, i_center]
        layers[2, mid_idx] = np.sqrt(r_center**2 - b_pc**2) * 2.0
        
        # 3. 前半球 (Front half)
        for i in range(mid_idx + 1, n_layer):
            i_in = i_center + i - (mid_idx + 1)
            i_out = i_in + 1
            r_out = shells[2, i_out]
            r_in = shells[2, i_in]
            layers[:, i] = shells[:, i_out]
            layers[2, i] = np.sqrt(r_out**2 - b_pc**2) - np.sqrt(r_in**2 - b_pc**2)
            
    return layers, n_layer

def calculate_central_density(mass_msun, alpha, r0_pc, rout_pc, cfg=default_config):
    """
    根据总质量和形状参数计算中心体密度 n_H2 (cm^-3)。
    
    参数:
    mass_msun : 拟合得到的总质量 [Msun]
    alpha     : 密度幂指数
    r0_pc     : 密度特征长度 [pc]
    rout_pc   : 云核总半径 [pc]
    """
    n_layers = 100 # 计算归一化常数时的解析精度
    rlayer_pc = rout_pc / n_layers
    r_pc = np.arange(n_layers) * rlayer_pc
    
    # 计算相对密度轮廓 (中心为 1.0)
    dens_rel = 1.0 / (1.0 + (r_pc / r0_pc)**alpha)
    
    # 计算每个壳层的物理体积 [cm^3]
    pc_to_cm = cfg.pc
    r_cm = r_pc * pc_to_cm
    dr_cm = rlayer_pc * pc_to_cm
    volumes = (4.0/3.0) * np.pi * ((r_cm + dr_cm)**3 - r_cm**3)
    
    # 计算中心密度为 1 cm^-3 时的总质量 [g]
    # 注意：mh2 是 H2 分子的质量，这里计算的是气体的总质量（已考虑气尘比转换等逻辑在模型内平衡）
    total_mass_rel_g = np.sum(dens_rel * volumes) * cfg.mh2
    
    # 中心密度 [cm^-3] = 真实总质量 [g] / 相对总质量 [g]
    nc = (mass_msun * cfg.msun) / total_mass_rel_g
    
    return nc

def calculate_mass_from_nc(nc_cm3, alpha, r0_pc, rout_pc, cfg=default_config):
    """
    根据中心密度和形状参数反求总质量 [Msun]。
    """
    n_layers = 100
    rlayer_pc = rout_pc / n_layers
    r_pc = np.arange(n_layers) * rlayer_pc
    
    dens_rel = 1.0 / (1.0 + (r_pc / r0_pc)**alpha)
    
    pc_to_cm = cfg.pc
    r_cm = r_pc * pc_to_cm
    dr_cm = rlayer_pc * pc_to_cm
    volumes = (4.0/3.0) * np.pi * ((r_cm + dr_cm)**3 - r_cm**3)
    
    total_mass_rel_g = np.sum(dens_rel * volumes) * cfg.mh2
    
    mass_g = nc_cm3 * total_mass_rel_g
    mass_msun = mass_g / cfg.msun
    
    return mass_msun
