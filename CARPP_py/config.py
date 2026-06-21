import numpy as np

class Config:
    """
    管理 CARPP 的全局物理常数和配置参数。
    替代 IDL 中的 common share, common parameter_range 和 setup_header.pro。
    """
    def __init__(self, grain_size_microns=0.1, distance_pc=414.0):
        # --- 基础物理常数 (CGS 单位制) ---
        self.h = 6.626e-27          # 普朗克常数 [erg s]
        self.k = 1.38e-16           # 玻尔兹曼常数 [erg/K]
        self.c = 3.0e10             # 光速 [cm/s]
        self.msun = 1.989e33        # 太阳质量 [g]
        self.jy = 1.0e-23           # 央斯基 [erg/s/cm^2/Hz]
        self.pc = 3.08568e18        # 秒差距 [cm]
        self.mh2 = 1.673e-24 * 2.0  # H2 分子质量 [g]
        self.pi = np.pi
        
        # --- 尘埃与环境参数 ---
        self.rho_d = 3.0            # 尘埃颗粒密度 [g/cm^3]
        self.q350 = 1.36e-4         # 350微米处的吸收效率 (Preibisch 1993)
        self.gdr = 100.0            # 气尘比 (Gas-to-Dust Ratio)
        
        # --- 派生常数 (基于输入参数) ---
        self.distance_pc = float(distance_pc)           # 距离 [pc]
        self.grain_size_cm = grain_size_microns * 1e-4  # 颗粒半径 [cm]
        self.distance_cm = float(self.distance_pc * self.pc) # 距离 [cm]
        self.arcsec_to_cm = self.distance_cm / 206265.0 # 1角秒对应的物理长度 [cm]

        # --- 拟合参数范围 (替代 setup_parameter_range.pro) ---
        self.t_index = 2.0          # 温度剖面指数 (Tindex)
        self.limit_range = False    # 是否限制参数范围
        self.delta_beta = 0.1
        self.delta_alpha = 0.2
        self.param_t_min = 3.0      # 最小允许温度 [K]
        self.delta_cloudmass = 1e9  # 云质量变化限制范围

    def update_derived(self, grain_size_microns, distance_pc):
        """当用户更改距离或颗粒大小时，更新派生常数"""
        self.distance_pc = float(distance_pc)
        self.grain_size_cm = grain_size_microns * 1e-4
        self.distance_cm = float(self.distance_pc * self.pc)
        self.arcsec_to_cm = self.distance_cm / 206265.0

# 实例化一个全局配置对象，方便其他模块调用
# 在 Python 中，这比 IDL 的 common 块更安全且易于调试
    @staticmethod
    def from_file(filepath):
        """从外部文本文件加载配置 (支持 Python 字面量语法)"""
        import os
        import ast
        settings = {}
        
        with open(filepath, 'r') as f:
            for line in f:
                # 1. 去除行内注释
                line = line.split('#', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                
                # 2. 分割键值
                key, value_str = line.split('=', 1)
                key = key.strip()
                value_str = value_str.strip()
                
                # 3. 使用 ast.literal_eval 还原 Python 字面量
                try:
                    val = ast.literal_eval(value_str)
                    # 如果是字符串，统一转小写清理空格
                    if isinstance(val, str):
                        val = val.strip().lower()
                    settings[key] = val
                except (ValueError, SyntaxError):
                    # 处理没加引号的字符串 (如 method = powell)
                    settings[key] = value_str.strip().lower()
        
        # 4. 路径处理：将相对路径转换为绝对路径 (基准是配置文件目录)
        base_dir = os.path.dirname(os.path.abspath(filepath))
        if 'image_files' in settings:
            files = settings['image_files']
            if isinstance(files, list):
                settings['image_files'] = [os.path.join(base_dir, f) for f in files]
            else:
                settings['image_files'] = [os.path.join(base_dir, files)]
        
        return settings

default_config = Config()
