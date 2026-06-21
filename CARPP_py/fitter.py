import numpy as np
from scipy.optimize import minimize, differential_evolution
from .config import default_config
from .models import get_core_model, convolve_gaussian

def _mcmc_log_prior(theta, bounds_arr):
    if np.all((theta >= bounds_arr[:, 0]) & (theta <= bounds_arr[:, 1])):
        return 0.0
    return -np.inf

def _mcmc_log_probability(theta, bounds_arr, fitter_instance):
    lp = _mcmc_log_prior(theta, bounds_arr)
    if not np.isfinite(lp):
        return -np.inf
    chi2 = fitter_instance.objective(theta)
    if chi2 >= 1e10:
        return -np.inf
    
    # 距离权重 weights_dist 被归一化到和为 1。因此 objective 返回的是每个像素的平均卡方误差。
    # 为了恢复统计学上真实的 log-likelihood 尺度，需要将卡方乘以当前的有效评估像素数 N_active。
    # 这样 log-likelihood 的绝对尺度才正确，MCMC walk step 才能得到正确的接受概率。
    n_active = np.sum(fitter_instance.weights_dist > 0)
    return lp - 0.5 * chi2 * n_active

def _grid_worker(args):
    """
    args = (fitter, chunk)
    """
    fitter, chunk = args
    results = []
    
    n_layers = 50
    rout_pc = fitter.rout_pc
    r_pc = (np.arange(n_layers) + 0.5) * (rout_pc / n_layers)
    pixel_size_pc = fitter.pixel_size_pc
    dist_pc = fitter.cfg.distance_pc
    pixel_size_arcsec = (pixel_size_pc / dist_pc) * 206265.0
    
    last_wave = fitter.wavelengths[-1]
    last_fwhm_pix = fitter.fwhm_beams[-1] / pixel_size_arcsec
    last_pad = int(np.ceil(last_fwhm_pix * 2.0))
    last_nx_model = fitter.nx + 2 * last_pad
    
    for beta, r0_layers, alpha, t0, t1, rt_layers in chunk:
        # Convert layers to pc
        r0_pc = r0_layers * (rout_pc / 50.0)
        rt_pc = rt_layers * (rout_pc / 50.0)
        
        dens_rel = 1.0 / (1.0 + (r_pc / r0_pc)**alpha)
        temp_prof = t1 + (t0 - t1) / (1.0 + (r_pc / rt_pc)**2.0)
        
        try:
            # 2. 用 M_ref = 1.0 生成最长波长处的 convolved model
            model_img_large, _ = get_core_model(
                temp_prof, dens_rel, 1.0, rout_pc, 
                fitter.cfg.grain_size_cm*1e4, beta, dist_pc,
                last_wave, n_layers, last_nx_model, pixel_size_pc,
                offset=fitter.offsets[-1], cfg=fitter.cfg
            )
            model_conv_large = convolve_gaussian(model_img_large, last_fwhm_pix)
            
            c_idx = (last_nx_model - 1) / 2.0
            icent = (fitter.nx - 1) / 2.0
            y_s = int(np.round(c_idx - icent))
            x_s = int(np.round(c_idx - icent))
            model_conv = model_conv_large[y_s:y_s+fitter.nx, x_s:x_s+fitter.nx]
            
            # Treating the longest-wavelength map as optically transparent:
            # cloud_mass = total_obs / total_model_1_0 (using weights_dist)
            total_obs = np.sum(fitter.obs_images[-1] * fitter.weights_dist)
            total_model = np.sum(model_conv * fitter.weights_dist)
            cloud_mass = total_obs / max(total_model, 1e-20)
            
            # Constrain cloud_mass
            if cloud_mass <= 0:
                cloud_mass = 1e-5
            elif cloud_mass > 1e10:
                cloud_mass = 1e10
                
            # 3. 计算所有波段的卡方值
            total_chi_sq = 0.0
            for i in range(len(fitter.wavelengths) - 1):
                fwhm_pix = fitter.fwhm_beams[i] / pixel_size_arcsec
                pad = int(np.ceil(fwhm_pix * 2.0))
                nx_model = fitter.nx + 2 * pad
                
                model_img_large_i, _ = get_core_model(
                    temp_prof, dens_rel, cloud_mass, rout_pc, 
                    fitter.cfg.grain_size_cm*1e4, beta, dist_pc,
                    fitter.wavelengths[i], n_layers, nx_model, pixel_size_pc,
                    offset=fitter.offsets[i], cfg=fitter.cfg
                )
                model_conv_large_i = convolve_gaussian(model_img_large_i, fwhm_pix)
                
                c_idx_i = (nx_model - 1) / 2.0
                icent_i = (fitter.nx - 1) / 2.0
                y_s_i = int(np.round(c_idx_i - icent_i))
                x_s_i = int(np.round(c_idx_i - icent_i))
                model_conv_i = model_conv_large_i[y_s_i:y_s_i+fitter.nx, x_s_i:x_s_i+fitter.nx]
                
                diff_sq_norm = (fitter.obs_images[i] - model_conv_i)**2 / fitter.variance_maps[i]
                band_chi = np.sum(diff_sq_norm * fitter.weights_dist)
                total_chi_sq += band_chi
                
            # 最后一波波段卡方并叠加
            diff_sq_norm_last = (fitter.obs_images[-1] - (model_conv * cloud_mass))**2 / fitter.variance_maps[-1]
            band_chi_last = np.sum(diff_sq_norm_last * fitter.weights_dist)
            total_chi_sq += band_chi_last
            
            # 整合积分 SED 的卡方贡献
            if fitter.sed_data is not None and len(fitter.sed_wavelengths) > 0:
                sed_chi_sum = 0.0
                if fitter.sed_width_arcsec is not None:
                    width_pc = (fitter.sed_width_arcsec / 206265.0) * dist_pc
                elif fitter.sed_width_pc is not None:
                    width_pc = fitter.sed_width_pc
                else:
                    width_pc = (89.0 / 206265.0) * dist_pc

                center_pix = (fitter.nx - 1) / 2.0
                y, x = np.indices((fitter.nx, fitter.nx))
                dx_pc = np.abs(x - center_pix) * pixel_size_pc
                dy_pc = np.abs(y - center_pix) * pixel_size_pc
                ap_mask = (dx_pc <= width_pc / 2.0) & (dy_pc <= width_pc / 2.0)
                if np.sum(ap_mask) == 0:
                    cx, cy = int(np.round(center_pix)), int(np.round(center_pix))
                    ap_mask[cy, cx] = True

                Dspec = 0.85
                for j in range(len(fitter.sed_wavelengths)):
                    sed_lam = fitter.sed_wavelengths[j]
                    model_sed_img, _ = get_core_model(
                        temp_prof, dens_rel, cloud_mass, rout_pc, 
                        fitter.cfg.grain_size_cm*1e4, beta, dist_pc,
                        sed_lam, n_layers, fitter.nx, pixel_size_pc,
                        offset=(0,0), cfg=fitter.cfg
                    )
                    fwhm_deg = (sed_lam * 1e-6 / Dspec) * (180.0 / np.pi)
                    fwhm_arcsec = fwhm_deg * 3600.0
                    fwhm_pix = fwhm_arcsec / pixel_size_arcsec
                    model_sed_conv = convolve_gaussian(model_sed_img, fwhm_pix)
                    
                    modflux = np.sum(model_sed_conv[ap_mask])
                    sed_chi_sum += (modflux - fitter.sed_fluxes[j])**2 / fitter.sed_variances[j]
                
                total_chi_sq += (sed_chi_sum / len(fitter.sed_wavelengths))

            results.append((total_chi_sq, [cloud_mass, beta, alpha, t0, t1, rt_pc, r0_pc]))
        except Exception as e:
            results.append((1e12, [0.0, beta, alpha, t0, t1, rt_pc, r0_pc]))
            
    return results

class CoreFitter:
    """
    CARPP 拟合引擎。
    """
    def __init__(self, image_data, wavelengths, fwhm_beams, subpixel_offsets=None, qz1=2.0, 
                 noise_rms=None, cal_error=0.1, rout_pc=0.2, pixel_size_pc=0.01, cfg=default_config,
                 sed_data=None, sed_noise_rms=0.0, sed_cal_error=0.1):
        self.obs_images = image_data  # [n_bands, nx, nx]
        self.wavelengths = wavelengths
        self.fwhm_beams = fwhm_beams
        self.offsets = subpixel_offsets if subpixel_offsets is not None else [(0,0)]*len(wavelengths)
        self.cfg = cfg
        self.nx = image_data.shape[1]
        self.qz1 = qz1
        self.rout_pc = rout_pc
        self.pixel_size_pc = pixel_size_pc
        self.weights_dist = self._calculate_weights(qz1=self.qz1)
        
        # 预计算误差方差图 (Variance Maps)
        # sigma^2 = rms^2 + (cal_error * I_obs)^2
        if noise_rms is None:
            self.noise_rms = np.ones(len(wavelengths)) * 0.1
        elif np.isscalar(noise_rms):
            self.noise_rms = np.ones(len(wavelengths)) * noise_rms
        else:
            noise_rms_arr = np.array(noise_rms)
            if noise_rms_arr.ndim == 0:
                self.noise_rms = np.ones(len(wavelengths)) * float(noise_rms_arr)
            elif len(noise_rms_arr) == 1:
                self.noise_rms = np.ones(len(wavelengths)) * noise_rms_arr[0]
            elif len(noise_rms_arr) != len(wavelengths):
                self.noise_rms = np.ones(len(wavelengths)) * noise_rms_arr[0]
            else:
                self.noise_rms = noise_rms_arr

        self.cal_error = cal_error
        self.variance_maps = []
        for i in range(len(wavelengths)):
            v_map = (self.noise_rms[i]**2) + (self.cal_error * self.obs_images[i])**2
            # 防止分母为 0
            v_map = np.maximum(v_map, 1e-10)
            self.variance_maps.append(v_map)

        # 辅助整合高积分 SED 信息进行联合拟合
        self.sed_data = sed_data
        self.sed_noise_rms = sed_noise_rms
        self.sed_cal_error = sed_cal_error
        if self.sed_data is not None:
            self.sed_wavelengths = np.array(sed_data['wavelengths'])
            self.sed_fluxes = np.array(sed_data['fluxes'])
            self.sed_width_arcsec = sed_data.get('width_arcsec')
            self.sed_width_pc = sed_data.get('width_pc')
            if 'errors' in sed_data and sed_data['errors'] is not None:
                self.sed_variances = np.array(sed_data['errors'])**2 + (self.sed_cal_error * self.sed_fluxes)**2
            else:
                self.sed_variances = (self.sed_noise_rms**2) + (self.sed_cal_error * self.sed_fluxes)**2
            self.sed_variances = np.maximum(self.sed_variances, 1e-10)

    def _calculate_weights(self, qz1=2.0):
        """计算距离权重矩阵 (Spatial Weighting)"""
        nx = self.nx
        center = (nx - 1) / 2.0
        y, x = np.indices((nx, nx))
        dist = np.sqrt((x - center)**2 + (y - center)**2)
        dist = np.maximum(dist, 1.0)
        
        w = 1.0 / (dist**qz1)
        w[dist > center] = 0.0
        # 归一化使得该波段的总权重和为 1，方便后续比较
        return w / (np.sum(w) + 1e-20)

    def objective(self, params, *args):
        # ... (参数约束逻辑保持不变)
        p = np.array(params)
        penalty = 0.0
        if p[0] <= 0: penalty += 1e10
        if p[3] <= 5: penalty += 1e10
        if p[4] <= 5: penalty += 1e10
        if p[5] <= 0.0001: penalty += 1e10
        if p[6] <= 0.0001: penalty += 1e10
        if penalty > 0: return penalty

        cloud_mass, beta, alpha, t0, t1, rt, r0 = p
        n_layers = 50
        rout_pc = self.rout_pc
        r_pc = (np.arange(n_layers) + 0.5) * (rout_pc / n_layers)
        dens_rel = 1.0 / (1.0 + (r_pc / r0)**alpha)
        temp_prof = t1 + (t0 - t1) / (1.0 + (r_pc / rt)**2.0)
        
        total_chi_sq = 0.0
        for i in range(len(self.wavelengths)):
            try:
                # 卷积到对应波段的分辨率
                # 使用从 Header 计算出的真实物理像素大小
                pixel_size_pc = self.pixel_size_pc
                dist_pc = self.cfg.distance_pc
                pixel_size_arcsec = (pixel_size_pc / dist_pc) * 206265.0
                fwhm_pix = self.fwhm_beams[i] / pixel_size_arcsec
                
                # 1. 确定计算网格大小 (nx_model)
                # 为防止卷积边际效应，我们需要给图像增加 Padding。
                # 同时，如果物理云核很大，也应该保证计算网格尽量覆盖云核。
                pad = int(np.ceil(fwhm_pix * 2.0))
                nx_model = self.nx + 2 * pad
                
                # 如果 rout_pc 对应的半径远大于当前网格，我们至少要在计算时包含它
                # 但为了计算量考虑，我们主要通过 pad 解决边际模糊问题
                # core_nx = int(np.ceil(2.0 * rout_pc / pixel_size_pc)) + 1
                # nx_model = max(nx_model, core_nx)
                
                model_img_large, _ = get_core_model(
                    temp_prof, dens_rel, cloud_mass, rout_pc, 
                    self.cfg.grain_size_cm*1e4, beta, self.cfg.distance_pc,
                    self.wavelengths[i], n_layers, nx_model, self.pixel_size_pc,
                    offset=self.offsets[i], cfg=self.cfg
                )
                
                # 2. 对大图进行卷积
                model_conv_large = convolve_gaussian(model_img_large, fwhm_pix)
                
                # 3. 裁剪回原始大小 self.nx x self.nx
                c_idx = (nx_model - 1) / 2.0
                icent = (self.nx - 1) / 2.0
                y_s = int(np.round(c_idx - icent))
                x_s = int(np.round(c_idx - icent))
                model_conv = model_conv_large[y_s:y_s+self.nx, x_s:x_s+self.nx]
                
                # 计算该像素的卡方项
                # (Obs - Model)^2 / Variance
                diff_sq_norm = (self.obs_images[i] - model_conv)**2 / self.variance_maps[i]
                
                # 结合距离权重进行求和
                # 每个波段贡献其平均加权误差，从而实现波段间的自动平衡
                band_chi = np.sum(diff_sq_norm * self.weights_dist)
                
                if np.isnan(band_chi): return 1e12
                total_chi_sq += band_chi
            except Exception as e:
                return 1e12
                
        # 2. 整合积分 SED 的卡方贡献（采用平衡权重的归一化卡方项）
        if self.sed_data is not None and len(self.sed_wavelengths) > 0:
            sed_chi_sum = 0.0
            
            # 模型中物理像素尺度与当前秒差距距离
            pixel_size_pc = self.pixel_size_pc
            dist_pc = self.cfg.distance_pc
            
            # 计算对应的物理光圈大小 (Square Aperture width in pc)
            if self.sed_width_arcsec is not None:
                width_pc = (self.sed_width_arcsec / 206265.0) * dist_pc
            elif self.sed_width_pc is not None:
                width_pc = self.sed_width_pc
            else:
                width_pc = (89.0 / 206265.0) * dist_pc

            # 预先生成积分的光圈掩膜 (Aperture Mask)
            center_pix = (self.nx - 1) / 2.0
            y, x = np.indices((self.nx, self.nx))
            dx_pc = np.abs(x - center_pix) * pixel_size_pc
            dy_pc = np.abs(y - center_pix) * pixel_size_pc
            ap_mask = (dx_pc <= width_pc / 2.0) & (dy_pc <= width_pc / 2.0)
            
            # 安全防线：如果光圈过小没有像素被选中，则保留最中心的一个像素
            if np.sum(ap_mask) == 0:
                cx, cy = int(np.round(center_pix)), int(np.round(center_pix))
                ap_mask[cy, cx] = True

            Dspec = 0.85 # 望远镜口径以米为单位
            for j in range(len(self.sed_wavelengths)):
                try:
                    sed_lam = self.sed_wavelengths[j]
                    
                    # 生成当前 SED 波段下的二维模型亮度图像，注意积分 SED 为中心对齐 (offset=0,0)
                    model_sed_img, _ = get_core_model(
                        temp_prof, dens_rel, cloud_mass, rout_pc, 
                        self.cfg.grain_size_cm*1e4, beta, dist_pc,
                        sed_lam, n_layers, self.nx, pixel_size_pc,
                        offset=(0,0), cfg=self.cfg
                    )
                    
                    # 对一维模型剖面按对应望远镜口径衍射限进行成像卷积
                    fwhm_deg = (sed_lam * 1e-6 / Dspec) * (180.0 / np.pi)
                    fwhm_arcsec = fwhm_deg * 3600.0
                    pixel_size_arcsec = (pixel_size_pc / dist_pc) * 206265.0
                    fwhm_pix = fwhm_arcsec / pixel_size_arcsec
                    
                    model_sed_conv = convolve_gaussian(model_sed_img, fwhm_pix)
                    
                    # 空间积分，求其射入观测光圈内的模型通量
                    modflux = np.sum(model_sed_conv[ap_mask])
                    
                    # 对该波长点贡献的卡方度量进行累加
                    term = (modflux - self.sed_fluxes[j])**2 / self.sed_variances[j]
                    sed_chi_sum += term
                except:
                    return 1e12
            
            # 将 SED 整体卡方除以其点数，将其转为归一化的一维光谱波段响应，
            # 类似上面每个二维频带各自归一化后的贡献。这使得 SED 能获得平衡的物理权重，
            # 既打破了简并，又使得二维图的对角和空隙特征依然得到良好的约束。
            total_chi_sq += (sed_chi_sum / len(self.sed_wavelengths))

        return total_chi_sq

    def fit(self, initial_guess=None, bounds=None, method='powell', workers=1, popsize=15, maxiter=100, tol=0.01, nwalkers=32, nsteps=1000):
        """
        执行拟合。
        method: 'powell', 'global' 或 'mcmc'
        popsize: 种群乘数 (仅用于 global)
        maxiter: 最大进化代数 (仅用于 global)
        tol: 收敛公差 (仅用于 global)
        nwalkers: MCMC 走步者数量 (仅用于 mcmc)
        nsteps: MCMC 迭代步数 (仅用于 mcmc)
        """
        if method == 'powell':
            print("Starting Powell minimization (Local Search)...")
            res = minimize(
                self.objective, 
                initial_guess, 
                method='Powell',
                bounds=bounds,
                options={'xtol': 1e-4, 'disp': True}
            )
            return res.x, res.fun, getattr(res, 'nit', 0)

        elif method == 'mcmc':
            import emcee
            import corner
            import matplotlib.pyplot as plt
            from multiprocessing import Pool
            import time

            print(f"Starting MCMC (Bayesian Sampling)...")
            print(f"Walkers: {nwalkers}, Steps: {nsteps}, Workers: {workers}")

            ndim = 7
            bounds_arr = np.array(bounds)

            # 1. 采用“联合两步法”方案：
            # 为保证绝对精度，省去可能陷入局部平庸势垒或低精度的自动化内部预寻优（如粗略 pre-fit）。
            # 最精细的科学工作流推荐：用户先通过高精度 'global' 模式寻找到极致精度全局最优点，
            # 随后将该最优参数填入 settings 的 'initial_guess' 中，运行 'mcmc' 模式。
            # MCMC 将直接以该精密最优点为中心，建立 Walkers 集群进行贝叶斯后验采样度量不确定性。
            print("MCMC initialized directly at the provided 'initial_guess' (highly recommended to be the output of a prior High-Precision GLOBAL run).")
            initial_p = np.array(initial_guess)

            # 2. 初始位置分布：根据各自物理参数边界的尺度 (bounds range) 设定抖动大小，比单纯固定 1e-4 更具物理关联
            bounds_range = bounds_arr[:, 1] - bounds_arr[:, 0]
            pos = initial_p + 1e-3 * bounds_range * np.random.randn(nwalkers, ndim)
            
            # 确保初始位置都在先验边界内
            for i in range(nwalkers):
                while not np.isfinite(_mcmc_log_prior(pos[i], bounds_arr)):
                    pos[i] = initial_p + 5e-4 * bounds_range * np.random.randn(ndim)

            pool = Pool(processes=workers if workers > 0 else None)
            sampler = emcee.EnsembleSampler(
                nwalkers, ndim, _mcmc_log_probability, 
                args=[bounds_arr, self], 
                pool=pool
            )
            
            start = time.time()
            sampler.run_mcmc(pos, nsteps, progress=True)
            end = time.time()
            pool.close()
            pool.join()

            print(f"MCMC finished in {end - start:.2f} seconds.")

            # 数据处理与可视化
            labels_raw = ["Mass", "Beta", "Alpha", "T0", "T1", "rt", "r0"]
            labels_nc = [r"$\beta$", r"$\rho_0$", r"$\alpha$", r"$r_0$", r"$T_0$", r"$T_1$", r"$r_{\rm t}$"]
            display_order_indices = [1, 0, 2, 6, 3, 4, 5]

            # 1. Trace Plot
            trace_samples = sampler.get_chain()
            fig, axes = plt.subplots(ndim, figsize=(10, 12), sharex=True)
            for i, orig_idx in enumerate(display_order_indices):
                ax = axes[i]
                ax.plot(trace_samples[:, :, orig_idx], "k", alpha=0.3)
                ax.set_ylabel(labels_nc[i] if i != 1 else "Mass [Msun]")
                ax.set_xlim(0, len(trace_samples))
                ax.yaxis.set_label_coords(-0.1, 0.5)
            axes[-1].set_xlabel("Step Number")
            plt.tight_layout()
            plt.savefig("mcmc_trace.png")
            print("Trace plot saved as 'mcmc_trace.png'.")

            # 2. Corner Plot (含 nc 转换)
            flat_samples = sampler.get_chain(discard=int(nsteps/2), thin=15, flat=True)
            from .physics import calculate_central_density
            
            final_samples = np.zeros_like(flat_samples)
            for i in range(len(flat_samples)):
                m, b, a, t0, t1, rt, r0 = flat_samples[i]
                nc = calculate_central_density(m, a, r0, self.rout_pc, cfg=self.cfg)
                final_samples[i] = [abs(b), nc, a, r0, t0, t1, rt]

            # Set font settings for larger text
            import matplotlib
            matplotlib.rcParams.update({
                'font.size': 14,
                'axes.labelsize': 16,
                'axes.titlesize': 14,
                'xtick.labelsize': 12,
                'ytick.labelsize': 12
            })

            fig = corner.corner(
                final_samples, labels=labels_nc, quantiles=[0.16, 0.5, 0.84],
                show_titles=True, title_fmt=".3f",
                label_kwargs={"fontsize": 16},
                title_kwargs={"fontsize": 14}
            )
            
            # 为外围坐标轴标签添加物理单位，而对角线处公式标题(Title)依然使用无单位的简洁形式
            axis_labels_with_units = [
                r"$\beta$",
                r"$\rho_0$ (cm$^{-3}$)",
                r"$\alpha$",
                r"$r_0$ (pc)",
                r"$T_0$ (K)",
                r"$T_1$ (K)",
                r"$r_{\rm t}$ (pc)"
            ]
            for ax in fig.axes:
                xlabel = ax.get_xlabel()
                ylabel = ax.get_ylabel()
                if xlabel in labels_nc:
                    idx = labels_nc.index(xlabel)
                    ax.set_xlabel(axis_labels_with_units[idx], fontsize=16)
                if ylabel in labels_nc:
                    idx = labels_nc.index(ylabel)
                    ax.set_ylabel(axis_labels_with_units[idx], fontsize=16)
                ax.tick_params(labelsize=12)
            plt.savefig("degeneracy_analysis.png")
            print("Corner plot saved as 'degeneracy_analysis.png'.")

            # 计算中位数作为结果
            mcmc_results = np.percentile(flat_samples, 50, axis=0)
            return mcmc_results, self.objective(mcmc_results), nsteps

        elif method == 'grid':
            import itertools
            from multiprocessing import Pool
            import time
            print("Starting Parallel Grid Search...")
            
            beta_grid = [-1.5, -1.6, -1.7, -1.8, -1.9, -2.0]
            r0_grid = [0.78, 1.56, 3.13, 6.25, 12.5, 25.0]
            alpha_grid = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
            t0_grid = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
            t1_grid = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
            rt_grid = [0.78, 1.56, 3.13, 6.25, 12.5, 25.0]
            
            all_combos = list(itertools.product(beta_grid, r0_grid, alpha_grid, t0_grid, t1_grid, rt_grid))
            n_combos = len(all_combos)
            print(f"Total grid combinations to evaluate: {n_combos}")
            
            # Chunking to balance multiprocessing overhead
            chunk_size = max(1, n_combos // (workers * 4 if workers > 0 else 16))
            chunks = [all_combos[i:i + chunk_size] for i in range(0, n_combos, chunk_size)]
            
            # Prepare pool args
            pool_args = [(self, chunk) for chunk in chunks]
            
            start = time.time()
            pool = Pool(processes=workers if workers > 0 else None)
            chunk_results = pool.map(_grid_worker, pool_args)
            pool.close()
            pool.join()
            end = time.time()
            print(f"Grid search completed in {end - start:.2f} seconds.")
            
            # Flatten results
            flat_results = []
            for res_list in chunk_results:
                flat_results.extend(res_list)
                
            # Find the best set of parameters
            best_chi2, best_p = min(flat_results, key=lambda x: x[0])
            
            return np.array(best_p), best_chi2, n_combos

        else:
            print(f"Starting Differential Evolution (Global Search)...")
            print(f"Workers: {workers}, Popsize multiplier: {popsize}, Max Generations: {maxiter}, Tol: {tol}")
            if bounds is None:
                raise ValueError("Global search requires 'bounds' parameter.")
            
            res = differential_evolution(
                self.objective,
                bounds=bounds,
                strategy='best1bin',
                maxiter=maxiter,
                popsize=popsize,
                tol=tol,
                mutation=(0.5, 1),
                recombination=0.7,
                disp=False,
                polish=True, # 开启最后的 L-BFGS-B 精细化
                updating='deferred' if workers != 1 else 'immediate',
                workers=workers
            )
            
        return res.x, res.fun, getattr(res, 'nit', 0)
