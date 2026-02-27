# Phase 2 — Problem Validation Report
# Research Type: C (Software Tool)
# Generated: Phase 2, Round 1 (convergence achieved)

═══════════════════════════════════════════════════════════════════
## Research Question (Validated)

**"HapGraph: 一个基于 F-统计量与群体间 IBD 共享统计量联合贝叶斯推断的
群体混合图（Population Admixture Graph）推断软件包，
同时估计混合比例（α）与混合时间（T，单位：代），
附 1kGP 2022 26 群体应用。"**

工具定位: 比 TreeMix 和 AdmixtureBayes 更先进 — 首个将 IBD 信息
纳入群体图推断的工具。

═══════════════════════════════════════════════════════════════════
## Round 1 Agent Verdicts

| Agent              | Verdict       | Core Concern |
|--------------------|---------------|-------------|
| Target User        | CONDITIONAL   | 运行时间需在数小时内；IBD 预处理需清晰文档；验证要扎实 |
| Journal Editor     | CONDITIONAL   | 需 ablation（F-only vs F+IBD）；timing 精度在模拟中验证；需第二个应用数据集 |
| Software Architect | CONDITIONAL   | 放弃全 RJMCMC，用 greedy 拓扑搜索 + 固定拓扑 MCMC；处理 F-stat/IBD 相关性 |

**Round 1 一致性评估:** 三个智能体条件高度一致，无根本性 FAIL → 收敛，无需 Round 2。

═══════════════════════════════════════════════════════════════════
## Key Insights from Deliberation

### 真实痛点 (Target User 确认)
- 现有工作流: TreeMix 推断拓扑 → ALDER/ROLLOFF 单独估计混合时间
- 这是两步走，不一致，难以做不确定性传播
- HapGraph 将两步合并为统一贝叶斯框架 → 真实痛点

### 新颖性 (Editor 确认)
- AdmixtureBayes (PLoS Genetics 2023) 是最近竞品：贝叶斯 + F-统计量，但 **无 IBD、无时间推断**
- HapGraph vs AdmixtureBayes 的差距就是 IBD 时间信息 → 核心创新点
- Reviewer 2 会问: "加了 IBD 真的提升了推断精度吗？" → 必须有 ablation

### 架构决策 (Software Architect 建议)
- **放弃全 RJMCMC** (N=26, K≤4 admixture 拓扑空间过大，单人6-12月不可行)
- **采用**: Greedy 拓扑搜索 (TreeMix-like) + PyMC NUTS 推断 α, T, N_e
- F-stat 与 IBD 统计量非独立，但互补：处理方式采用块对角协方差矩阵
- N_e 可固定（来自 PSMC 外部估计）以减少参数数量

═══════════════════════════════════════════════════════════════════
## Refined Architecture (Post-Phase-2)

```
INPUT: 相位化 VCF + 群体标签
       ↓
STEP 1: PREPROCESSING
  - ADMIXTOOLS2/admixr: 计算 F2, F3, F4 统计量
  - hap-ibd:           计算群体间 IBD 段 (≥ 2cM)
  - Custom script:     聚合 IBD → L̄, r per 群体对
       ↓
STEP 2: TOPOLOGY SEARCH (greedy)
  - 从 NJ/Neighbor-Joining 树出发
  - 贪婪添加 admixture 边 (逐步最大化联合似然)
  - 输出 top-5 候选拓扑
       ↓
STEP 3: BAYESIAN INFERENCE (PyMC, NUTS)
  - 固定拓扑，推断 α, T, [N_e]
  - 联合似然: L(F-stats, IBD | α, T, N_e)
  - 块对角协方差 (F-stat 块 + IBD 块)
  - 3 链并行, R̂ < 1.05 收敛标准
       ↓
STEP 4: MODEL COMPARISON
  - BIC 跨候选拓扑
  - 输出最优拓扑 + 参数后验
       ↓
OUTPUT: 
  - 最优混合图 (可视化, Newick-like 格式)
  - α 和 T 的后验分布 + 95% CI
  - 诊断报告 (trace plots, ESS, R̂)
```

═══════════════════════════════════════════════════════════════════
## Required Deliverables (for paper acceptance)

| # | Deliverable | Why Required |
|---|------------|--------------|
| 1 | msprime 模拟基准 (already-known topology+timing) | 验证方法准确性 |
| 2 | **Ablation**: F-only vs F+IBD 精度对比 | Reviewer 2 核心关切 |
| 3 | AdmixtureBayes 比较 (相同数据) | 同类最强竞品 |
| 4 | TreeMix 比较 (相同数据) | 领域标准基线 |
| 5 | 1kGP 2022 26群体应用 | 真实数据展示 |
| 6 | 第二应用数据集 (HGDP 或模拟已知真值场景) | 编辑要求 |
| 7 | 软件包 (pip/conda 安装, 文档, 示例) | 工具论文必须 |

═══════════════════════════════════════════════════════════════════
## Conditions Accepted (CONDITIONAL → treated as PASS)

- [ ] **COND-1**: 拓扑搜索用 greedy，不用全 RJMCMC  ← 在 Phase 3 方法设计中锁定
- [ ] **COND-2**: 提供 F-only ablation 实验            ← 在 Phase 3 实验计划中锁定
- [ ] **COND-3**: 模拟数据 timing 精度验证              ← Phase 3 评估协议中锁定
- [ ] **COND-4**: 处理 IBD 与 F-stat 相关性 (块对角)   ← Phase 3 方法设计中锁定
- [ ] **COND-5**: 添加第二应用数据集                    ← Phase 3 实验计划中锁定

═══════════════════════════════════════════════════════════════════
## Project Intent Classification (Step 2)

Primary Intent: **Engineering** (build a practical, deployable tool)
Secondary:      **Quality** (improve reliability via uncertainty quantification)

Research Type: C (confirmed by user)
