dsRNA-Forge Windows 傻瓜式使用说明

启动：
1. 解压整个发布包。
2. 双击 dsRNA-Forge.exe。
3. 不需要安装 Python，不需要运行 pip，不需要配置环境变量。

第一次试运行：
1. 点击 Browse...。
2. 选择 demo_data\test_transcriptome.fa。
3. 点击 Load & Index。
4. 在 Target source 里选择目标来源：
   - Search in transcriptome：在转录组里搜索并选择。
   - Paste target sequence：直接粘贴一条 target 序列。
   - Upload target FASTA：上传单条 target FASTA。
5. 选择 siRNA、DsiRNA、Long dsRNA 或 sgRNA for SpCas9。
6. 在 Design confidence 里选择：
   - Strict：脱靶更严格，结果更少。
   - Balanced：默认推荐。
   - Relaxed：目标难设计或 0 结果时使用。
7. 可选：点击 Add Background FASTA... 加入宿主、近缘非靶标、益虫等额外转录组作为脱靶背景。
8. 普通用户不需要展开 Show advanced settings，长度、GC、规则和 CPU 核心数会自动使用默认值。
9. 点击 Start Design。
10. 在 Results 页面查看结果；点选某一行可看推荐理由、脱靶验证片段、区域图和引物/寡核苷酸。
11. 使用 Export CSV / Export Excel / Export FASTA / Export Report / Export Primers 导出。

再次使用同一个转录组：
1. 打开软件后，在 Saved transcriptomes 下拉框选择之前加载过的转录组。
2. 点击 Load Saved。
3. 不需要重新 Browse；软件会从本地缓存载入转录组和脱靶索引。
4. 点击 Manage Cache 可以重命名、删除或清空已入库转录组。

输入文件：
- 支持 .fa / .fasta / .fna。
- FASTA 标题行以 > 开头。
- DNA 序列中的 T 会自动转成 RNA 的 U。

说明：
- 当前包内置 PyQt6、Biopython、pandas、openpyxl 和 ViennaRNA。
- ViennaRNA 2.7.2 已通过运行时自检。
- 最终 Windows 包必须通过真实 RNAup CLI 自检；运行 dsRNA-Forge.exe --check-runtime 应显示 RNAup method: RNAup-cli。
- 如果结果为 0，软件会弹出可能原因和建议，例如目标太短、GC 过滤太窄或转录组里存在高度相似序列。
- 结果表会显示 Off-target Risk、Risk Score、Top Risk Targets 和 Validation Direction，用于安排后续 BLAST/qPCR/片段比对验证。
- Export CSV / Excel / Report 会明确列出 Top Off-target Genes，也就是潜在脱靶基因/转录本、风险分和主要风险原因。
- 软件会为 Top 候选尝试 RNAup 精筛；如果开发环境没有 RNAup CLI，会明确记录为 RNAduplex fallback，不会冒充 RNAup 精确结果。正式 Windows 包不应出现这个 fallback。
- 默认结果是非冗余推荐：相邻 1 bp 滑窗、重叠度很高的候选会折叠成一条代表结果。
- Cluster Size 表示这一条代表结果下面折叠了多少个相似候选；数值越大，说明附近有更多等价替代窗口。
- sgRNA 模式按 SpCas9 20nt spacer + NGG PAM 扫描正负链，结果会显示 PAM、strand、cut site、Top off-target 和验证方向。
- sgRNA 模式会导出 BbsI 兼容克隆寡核苷酸，以及切点附近 genotyping PCR 引物；高风险 off-target 建议做靶向扩增测序/ICE/TIDE/Sanger。
- Export Report 会生成实验验证报告，包含推荐理由、Top 脱靶验证对象、匹配片段、T7 引物、sgRNA oligo/分型引物和方法说明。
- Export Primers 会导出 T7 引物、sgRNA 克隆寡核苷酸和分型 PCR 引物订购 CSV。
- Save Project / Open Project 支持 .dsforge_project 项目文件，保存 target、参数和当前结果。

如果双击没反应：
1. 确认不是在 zip 压缩包内部直接运行。
2. 先完整解压到桌面。
3. 再双击 dsRNA-Forge.exe。
