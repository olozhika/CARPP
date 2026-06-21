import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
from multiprocessing import Pool
import CARPP_py
from CARPP_py.preprocessing import preprocess_images
from CARPP_py.fitter import CoreFitter
from CARPP_py.config import Config
import time

def log_prior(theta, bounds):
    if np.all((theta >= bounds[:, 0]) & (theta <= bounds[:, 1])):
        return 0.0
    return -np.inf

def log_likelihood(theta, fitter):
    # CARPP 的 objective 返回的是 chi_square
    # 似然为 exp(-0.5 * chi2)
    chi2 = fitter.objective(theta)
    if chi2 >= 1e10: # 处理 penalty
        return -np.inf
    return -0.5 * chi2

def log_probability(theta, fitter, bounds):
    lp = log_prior(theta, bounds)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta, fitter)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll

def run_mcmc_analysis(config_path, nwalkers=32, nsteps=2000, nthreads=None):
    """
    使用 MCMC 进行参数简并性分析
    """
    # 1. 加载设置与预处理 (复用 CARPP 逻辑)
    print(f"Loading configuration from: {config_path}")
    settings = Config.from_file(config_path)
    
    # 初始化物理环境配置
    cfg = Config(
        grain_size_microns=settings.get('grain_radius_microns', 0.1),
        distance_pc=settings.get('distance_pc', 414.0)
    )

    # 获取预处理所需的关键参数
    n_cells = int(settings['n_cells'])
    default_center = [(n_cells - 1) / 2.0, (n_cells - 1) / 2.0]
    
    # 正确调用 preprocess_images
    obs_images, subpixel_offsets, pixel_size_arcsec = preprocess_images(
        settings['image_files'], 
        settings['fwhm_arcsec'], 
        settings['back_levels'], 
        n_cells,
        center_type=settings.get('center_type', 'pixel'),
        center_val=settings.get('center_val', default_center),
        cfg=cfg
    )
    
    # 2. 拟合引擎准备
    qz1 = settings.get('qz1', 2.0)
    noise_rms = settings.get('noise_rms', [0.1]*len(settings['image_files']))
    cal_error = settings.get('cal_error', 0.1)
    
    distance_pc = cfg.distance_pc
    pixel_size_pc = (pixel_size_arcsec / 206265.0) * distance_pc
    
    # 智能处理 rout_pc
    rout_pc = settings.get('rout_pc')
    if rout_pc is None or str(rout_pc).lower() == 'none':
        rout_pc = ((n_cells - 1) / 2.0) * pixel_size_pc

    # 实例化 Fitter
    fitter = CoreFitter(
        obs_images, 
        settings['wavelengths_microns'], 
        settings['fwhm_arcsec'], 
        subpixel_offsets=subpixel_offsets, 
        qz1=qz1, 
        noise_rms=noise_rms, 
        cal_error=cal_error, 
        rout_pc=rout_pc, 
        pixel_size_pc=pixel_size_pc, 
        cfg=cfg
    )

    # 3. 参数与状态初始化
    # 参数顺序: [Mass, Beta, Alpha, T0, T1, rt, r0]
    bounds = np.array(settings['bounds'])
    ndim = 7
    # 以初始猜测为中心进行微小抖动初始化
    initial_p = np.array(settings.get('initial_guess', [30.0, -1.6, 2.0, 30.0, 15.0, 0.05, 0.05]))
    pos = initial_p + 1e-4 * np.random.randn(nwalkers, ndim)
    
    # 确保初始位置在边界内
    for i in range(nwalkers):
        while not np.isfinite(log_prior(pos[i], bounds)):
            pos[i] = initial_p + 1e-3 * np.random.randn(ndim)

    # 4. 开始采样 (并行)
    print(f"Starting MCMC with {nwalkers} walkers, {nsteps} steps...")
    
    # 诊断性检查：计算初始位置和边界附近的 chi2 差异
    chi2_initial = fitter.objective(initial_p)
    # 稍微挪动一下参数看 chi2 变化
    p_test = initial_p.copy()
    p_test[0] *= 1.1 # 增加 10% 质量
    chi2_test = fitter.objective(p_test)
    
    print(f"Initial Chi2: {chi2_initial:.2f}")
    print(f"Chi2 after 10% mass change: {chi2_test:.2f}")
    print(f"Delta Chi2: {abs(chi2_test - chi2_initial):.2f}")
    
    if abs(chi2_test - chi2_initial) < 1.0:
        print("WARNING: The likelihood surface is extremely flat. MCMC will perform a random walk.")
        print("ADVICE: Decrease 'noise_rms' or 'cal_error' in your settings file to make the fit more sensitive.")

    with Pool(processes=nthreads) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool, args=(fitter, bounds))
        start = time.time()
        sampler.run_mcmc(pos, nsteps, progress=True)
        end = time.time()
        print(f"MCMC finished in {end - start:.2f} seconds.")

    # 6. 数据导出与可视化
    # 目标顺序: ["Beta", "n_c [cm^-3]", "Alpha", "r0", "T0", "T1", "rt"]
    labels = ["Beta", "n_c [cm^-3]", "Alpha", "r0", "T0", "T1", "rt"]
    
    # 6.1 生成 Trace Plot (检查收敛性)
    # 注意：Trace 图为了效率展示原始 Mass，对应目标顺序中的 n_c 位置
    trace_samples = sampler.get_chain()
    fig, axes = plt.subplots(ndim, figsize=(10, 12), sharex=True)
    
    # 原始索引与展示顺序的对应关系
    # 原始: [0:Mass, 1:Beta, 2:Alpha, 3:T0, 4:T1, 5:rt, 6:r0]
    # 目标: [1:Beta, 0:Mass, 2:Alpha, 6:r0, 3:T0, 4:T1, 5:rt]
    display_order_indices = [1, 0, 2, 6, 3, 4, 5]
    
    for i, orig_idx in enumerate(display_order_indices):
        ax = axes[i]
        ax.plot(trace_samples[:, :, orig_idx], "k", alpha=0.3)
        ax.set_ylabel(labels[i] if i != 1 else "Mass [Msun]")
        ax.set_xlim(0, len(trace_samples))
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("Step Number")
    plt.tight_layout()
    plt.savefig("mcmc_trace.png")
    print("Trace plot saved as 'mcmc_trace.png' (Note: nc row displays Mass for speed).")

    # 6.2 生成 Corner Plot
    # 丢弃前 1/2 的阶段作为 Burn-in
    flat_samples = sampler.get_chain(discard=int(nsteps/2), thin=15, flat=True)
    
    # 参数转换与重排
    from CARPP_py.physics import calculate_central_density
    print("Transforming parameters and reordering for Corner plot...")
    
    final_samples = np.zeros_like(flat_samples)
    for i in range(len(flat_samples)):
        m, b, a, t0, t1, rt, r0 = flat_samples[i]
        nc = calculate_central_density(m, a, r0, rout_pc, cfg=cfg)
        # 顺序: ["Beta", "n_c [cm^-3]", "Alpha", "r0", "T0", "T1", "rt"]
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
        final_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
        show_titles=True, title_fmt=".3f",
        label_kwargs={"fontsize": 16},
        title_kwargs={"fontsize": 14}
    )
    
    # Increase tick labels size for all axes
    for ax in fig.axes:
        ax.tick_params(labelsize=12)
        # Check and update labels to ensure they are large
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=16)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=16)

    plt.savefig("degeneracy_analysis.png")
    print("Corner plot saved as 'degeneracy_analysis.png'")
    
    # 打印参数不确定度
    for i in range(ndim):
        mcmc = np.percentile(final_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(f"{labels[i]}: {mcmc[1]:.4f} (+{q[1]:.4f}, -{q[0]:.4f})")

    return final_samples

if __name__ == "__main__":
    # 调用示例
    path = 'data/carpp_settings.txt'
    samples = run_mcmc_analysis(path, nwalkers=32, nsteps=1000)