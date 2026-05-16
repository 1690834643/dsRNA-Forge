# dsRNA-Forge

**昆虫/植物 dsRNA/siRNA/DsiRNA 与 SpCas9 sgRNA 设计工具**

离线可用的 Windows 桌面应用，支持自定义转录组、转录组内搜索目标、自定义 target 序列、ViennaRNA RNAduplex/RNAcofold/RNAup CLI 热力学评估、SpCas9 sgRNA 设计和多核并行。

> 本仓库保存源码、测试、示例数据和打包配置；Windows 便携 exe 属于构建产物，不直接提交到仓库。实验前请把 Top off-target 结果作为验证方向，而不是“绝对无脱靶”的承诺。

---

## 核心特性

- **完全自定义转录组**：支持任意物种 FASTA 上传，不依赖内置数据库
- **目标序列三种来源**：可在转录组里搜索选择，也可直接粘贴 target 序列或上传单条 target FASTA
- **转录组入库缓存**：第一次加载后写入本地缓存，后续可直接从 Saved transcriptomes 选择
- **离线开箱即用**：Windows 双击运行，无需安装 Python/ViennaRNA/BLAST
- **长 dsRNA 专用设计逻辑**：预测 Dicer 切割位点，评估 siRNA pool 整体质量
- **热力学脱靶评估**：内嵌 ViennaRNA；Windows 包自检要求真实 `RNAup-cli` 可用，Top 候选会记录 RNAup 方法
- **sgRNA 设计模块**：支持 SpCas9 20nt spacer + NGG PAM，扫描正负链，输出 PAM、strand、cut site、on-target score
- **sgRNA 脱靶排序**：对 NGG-adjacent 近似位点做 mismatch 风险排序，给出 Top off-target 和扩增测序验证方向
- **sgRNA 实验寡核苷酸**：导出 pX330/lentiCRISPR 常见 BbsI 克隆寡核苷酸，以及切点附近基因分型 PCR 引物
- **无结果诊断**：当候选为 0 或全部被过滤时，会提示可能原因和下一步调整建议
- **脱靶风险排序**：显示风险等级、风险分、Top 风险转录本和验证方向，导出文件同步包含这些列
- **非冗余推荐**：默认折叠相邻 1 bp 滑窗候选，展示独立代表结果和 Cluster Size
- **实验决策面板**：点选结果即可查看推荐理由、扣分原因、脱靶验证片段、区域图和后续验证建议
- **引物/寡核苷酸设计**：Long dsRNA 自动生成普通 PCR primer 和 T7 promoter primer；sgRNA 自动生成克隆 oligo 和分型 PCR primer
- **多背景脱靶**：支持额外加载宿主、近缘非靶标或益虫转录组作为脱靶背景
- **项目文件和缓存管理**：支持 `.dsforge_project` 保存/打开，以及已入库转录组重命名、删除和清理
- **三档傻瓜模式**：Strict / Balanced / Relaxed，普通用户无需理解规则和核心数
- **多核并行**：充分利用现代 CPU，批量设计不卡顿
- **四种设计模式**：
  - siRNA 模式（21nt）
  - DsiRNA 模式（27nt）
  - 长 dsRNA 模式（200-1000bp）
  - sgRNA 模式（SpCas9, 20nt + NGG）

## 技术栈

| 组件 | 选型 | 版本 |
|------|------|------|
| 语言 | Python | 3.11+ |
| GUI | PyQt6 | 6.5+ |
| 序列处理 | Biopython | 1.81+ |
| 数值计算 | numpy, pandas | 最新稳定版 |
| 热力学 | ViennaRNA Python API + RNAup CLI | 2.7.2 |
| 数据存储 | SQLite | 标准库 |
| 打包 | PyInstaller | 6.0+ |

## 项目结构

```
dsRNA_Forge/
├── main.py                      # 程序入口
├── config.json                  # 默认配置
├── requirements.txt             # 依赖清单
├── dsRNA-Forge.spec             # PyInstaller 打包配置
├── dsforge/                     # 主包
│   ├── core/                    # 计算引擎（零 Qt 依赖）
│   │   ├── sequence.py          # FASTA 加载、索引、候选生成
│   │   ├── scoring/             # 多规则评分引擎
│   │   │   ├── base.py          # ScoringRule ABC + 注册表
│   │   │   ├── reynolds.py      # Reynolds (2004)
│   │   │   ├── ui_tei.py        # Ui-Tei (2004)
│   │   │   ├── amarzguioui.py   # Amarzguioui (2004)
│   │   │   ├── hsieh.py         # Hsieh (2004)
│   │   │   ├── jagla.py         # Jagla (2005)
│   │   │   └── consensus.py     # 共识评分器 (PMC5357899, 2017)
│   │   ├── dicer.py             # Dicer 切割位点预测
│   │   ├── offtarget.py         # dsRNA/siRNA 脱靶筛查
│   │   ├── sgrna.py             # SpCas9 sgRNA 扫描、脱靶排序、寡核苷酸
│   │   └── thermodynamics.py    # ViennaRNA 封装
│   ├── controller/              # 控制/服务层
│   │   ├── design_task.py       # 设计任务（单线程）
│   │   ├── design_task_parallel.py  # 设计任务（多进程）
│   │   ├── scheduler.py         # 多进程调度器
│   │   └── exporter.py          # 结果导出
│   ├── database/                # 数据层
│   │   ├── schema.py            # SQLite Schema
│   │   └── manager.py           # CRUD 封装
│   ├── gui/                     # GUI 层（PyQt6）
│   │   ├── main_window.py
│   │   ├── transcript_panel.py
│   │   ├── config_panel.py
│   │   ├── progress_panel.py
│   │   ├── results_panel.py
│   │   ├── history_panel.py
│   │   └── workers.py
│   └── utils/
│       └── vienna_loader.py     # 跨平台 ViennaRNA 加载
└── tests/                       # 测试
    ├── test_cli.py              # CLI 核心引擎测试
    ├── test_parallel.py         # 多进程测试
    ├── test_integration.py      # 端到端集成测试
    └── test_gui_headless.py     # GUI 组件测试
```

## 快速开始

### Windows 用户

如果你只想直接使用软件，请下载项目发布页中的 Windows portable 包，解压后双击 `启动 dsRNA-Forge.bat` 或 `dsRNA-Forge.exe`。源码仓库本身不包含 exe。

### 开发环境（Linux）

```bash
# 1. 创建虚拟环境
python3.11 -m venv venv
source venv/bin/activate

# 2. 安装依赖（推荐使用清华镜像）
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple

# 3. 验证安装
python -c "import RNA; print('ViennaRNA:', RNA.__version__)"
python -c "from PyQt6.QtWidgets import QApplication; print('PyQt6 OK')"

# 4. 运行 CLI 测试
python test_cli.py
python test_parallel.py
python test_integration.py

# 5. 启动 GUI
python main.py
```

### Windows 打包

```bash
# 1. 在 Windows 构建机上安装官方 ViennaRNA Windows 包
#    下载地址: https://www.tbi.univie.ac.at/RNA/

# 2. 设置环境变量
set VIENNA_DLL_DIR=C:\Program Files\ViennaRNA

# 3. 打包
pyinstaller dsRNA-Forge.spec

# 4. 输出目录
#    dist/dsRNA-Forge/
#    ├── dsRNA-Forge.exe
#    └── _internal/
```

## 使用说明

1. **加载转录组**：点击 "Browse..." 选择 FASTA 文件，点击 "Load & Index"
2. **选择目标序列**：
   - "Search in transcriptome"：在已加载转录组中搜索并选择目标 ID
   - "Paste target sequence"：直接粘贴一条 target 序列；DNA/RNA 均可
   - "Upload target FASTA"：上传单条 target FASTA
3. **配置参数**：
   - 普通用户只需要选择设计模式（siRNA / DsiRNA / Long dsRNA / sgRNA）和 Strict / Balanced / Relaxed
   - 长度范围、GC、评分规则和 CPU 核心数放在 "Show advanced settings" 中，默认隐藏
4. **额外背景**：可选加载宿主、近缘非靶标或益虫 FASTA 作为 off-target 背景
5. **开始设计**：点击 "Start Design"
6. **查看结果**：在 Results 标签页查看排序后的候选列表，点选一行查看解释和验证片段
7. **导出**：点击 Export CSV / Excel / FASTA / Report / Primers 保存结果

## 设计规则说明

| 规则 | 来源 | 标准数 | 特点 |
|------|------|--------|------|
| Reynolds | Dharmacon, 2004 | 8 | 经验规则，cutoff ≥ 6 |
| Ui-Tei | Ui-Tei et al., 2004 | 4 |  antisense 链特异性 |
| Amarzguioui | Amarzguioui et al., 2004 | 3 | 位置偏好 |
| Hsieh | Hsieh et al., 2004 | 3 | 位置特征 + 热力学 |
| Jagla | Jagla et al., 2005 | 4套 | 按 GC 含量分决策树 |
| **Consensus** | PMC5357899, 2017 | 12+ | **默认引擎**，多规则整合 |

## ViennaRNA 热力学说明

Windows 发布包随 exe 打入官方 ViennaRNA CLI，运行时自检要求 `RNAup.exe` 可被自动发现并返回 `RNAup-cli`。源码开发环境如果没有安装 CLI，测试仍允许明确标记的 fallback，但最终 Windows 包不能用 fallback 冒充 RNAup 精确结果。

| 工具 | 当前用途 | 说明 |
|------|----------|------|
| RNAduplex | 候选结合能、种子区热力学脱靶预筛 | 已启用 |
| RNAcofold | 两分子联合结构/能量计算接口 | 已封装，可用于后续更精细评估 |
| RNAup | 含 mRNA 结构打开惩罚的更精确模式 | Windows 发布包强制自检 `RNAup-cli`；导出列显示实际方法 |

导出列中 `rnaup_method=RNAup-cli` 的结果代表实际 RNAup CLI 精筛；若在开发环境中看到 fallback 标记，表示该环境没有可用的 RNAup CLI，不应作为最终交付包使用。

## sgRNA 模块说明

- 默认标准：SpCas9，20nt spacer，NGG PAM。
- 扫描范围：目标序列正链 NGG 和反链 CCN 都会扫描。
- 排序依据：on-target heuristic、GC 区间、U6 终止信号、PAM-proximal seed 区域和 off-target 风险。
- 脱靶输出：Top off-target、mismatch 数、PAM、strand、risk score 和验证方向。
- 实验输出：BbsI 兼容 forward/reverse cloning oligo，以及切点附近 genotyping PCR primer。
- 注意：当前离线模型是可解释 heuristic，不宣称等同 CRISPOR/CRISPRdirect/CHOPCHOP 的完整在线数据库评分；高价值实验仍建议对 Top off-target 做靶向扩增测序。

## 实验交付输出

- `Export Report`：多 sheet Excel，包含 Recommendations、OffTarget Validation、Primers 和 Methods；dsRNA/siRNA/DsiRNA 会显式列出 Top Off-target Genes、风险原因、匹配片段和建议验证动作。
- `Export Primers`：Long dsRNA 导出 T7 引物；sgRNA 导出克隆 oligo 和分型 PCR primer。
- `Save Project` / `Open Project`：保存或恢复 `.dsforge_project`，包含当前 target、参数、缓存引用和结果。

## 仓库内容

- `dsforge/`：核心算法、GUI、数据库和导出逻辑
- `demo_data/`：小型 FASTA 示例
- `test_*.py`：交付级回归测试和集成测试
- `dsRNA-Forge.spec`：PyInstaller Windows 打包配置
- `WINDOWS_USER_GUIDE.md` / `README_使用说明.txt`：面向 Windows 用户的操作说明

## 测试覆盖

```
Test 1: Sequence Utilities          ✓
Test 2: Scoring Engine (6 Rules)    ✓
Test 3: Dicer Cleavage Prediction   ✓
Test 4: Database                    ✓
Test 5: siRNA Mode                  ✓
Test 6: Long dsRNA Mode             ✓
Test 7: DsiRNA Mode                 ✓
Test 8: ViennaRNA Thermodynamics    ✓
Test 9: Result Exporter             ✓
Test 10: SpCas9 sgRNA Design        ✓
Parallel siRNA (4 cores)            ✓
Parallel Long dsRNA (4 cores)       ✓
GUI Components (Headless)           ✓
End-to-End Integration              ✓
```

## 许可证

MIT License

## 参考文献

1. Reynolds et al. (2004) — Rational siRNA design for RNA interference
2. Ui-Tei et al. (2004) — Guidelines for the selection of highly effective siRNA sequences
3. Jagla et al. (2005) — Sequence characteristics of functional siRNAs
4. PMC5357899 (2017) — Predicting siRNA efficacy based on multiple selective criteria
5. ViennaRNA Package 2.0 — https://www.tbi.univie.ac.at/RNA/
6. CRISPOR — Haeussler et al. (2016), Genome Biology
7. CRISPRdirect — Naito et al. (2015), Bioinformatics
8. CHOPCHOP v2 — Labun et al. (2016), Nucleic Acids Research
9. Doench et al. (2016), Nature Biotechnology
