# dsRNA-Forge Windows 使用说明

## 直接启动
1. 解压发布包。
2. 双击 `dsRNA-Forge.exe`。
3. 点击 `Browse...` 选择自己的转录组 FASTA 文件，或先用 `demo_data/test_transcriptome.fa` 试运行。
4. 点击 `Load & Index`。
5. 在 `Target source` 中选择目标来源：
   - `Search in transcriptome`：在转录组里搜索并选择 target。
   - `Paste target sequence`：直接粘贴一条 target 序列，DNA/RNA 均可。
   - `Upload target FASTA`：上传单条 target FASTA。
6. 选择设计模式：`siRNA`、`DsiRNA`、`Long dsRNA` 或 `sgRNA for SpCas9`。
7. 选择 `Design confidence`：
   - `Strict`：优先低脱靶风险，结果更少。
   - `Balanced`：默认推荐。
   - `Relaxed`：目标难设计或无结果时使用。
8. 可选：用 `Add Background FASTA...` 加入宿主、近缘非靶标、益虫等额外脱靶背景。
9. 普通用户不用展开 `Show advanced settings`；长度、GC、规则和 CPU 核心数已经有默认值。
10. 点击 `Start Design`。
11. 在右侧结果表查看候选；点选一行可查看解释、脱靶验证片段、区域图和引物/寡核苷酸。
12. 用 `Export CSV`、`Export Excel`、`Export FASTA`、`Export Report` 或 `Export Primers` 导出。

## 再次使用同一个转录组
- 第一次 `Load & Index` 后，软件会把转录组写入本地缓存和已保存列表。
- 下次打开软件，在 `Saved transcriptomes` 里选中之前的转录组，点击 `Load Saved` 即可。
- 同一转录组的脱靶 k-mer 索引也会缓存；后续预测不同基因会复用该索引。
- `Manage Cache` 可重命名、删除或清理已入库转录组。

## 输入文件要求
- 支持 `.fa`、`.fasta`、`.fna`。
- FASTA 标题行以 `>` 开头。
- 序列可以是 DNA 或 RNA；软件会自动把 `T` 转换为 `U`。

## 不需要用户安装的内容
- 不需要安装 Python。
- 不需要运行 `pip install`。
- 不需要配置环境变量。
- 不需要单独安装 PyQt6、Biopython、pandas 或 ViennaRNA。

## 当前热力学说明
- 发布包内置 ViennaRNA Python 绑定和官方 ViennaRNA CLI。
- `dsRNA-Forge.exe --check-runtime` 必须显示 `RNAup method: RNAup-cli`，否则这个 Windows 包不合格。
- Top 候选会尝试 RNAup 精筛；导出文件中的 `RNAup Method` 显示实际使用的方法。
- 源码开发环境缺少 RNAup CLI 时可能出现 RNAduplex fallback 标记，但最终 Windows 便携包不能把 fallback 当成 RNAup 精确结果交付。

## 脱靶风险排序
- `Risk` / `Risk Score`：按 20bp 连续匹配、16bp 连续匹配、7nt seed 命中等信息汇总。
- `Top Risk Targets`：最需要优先验证的潜在脱靶转录本。
- `Top Off-target Genes`：导出报告里更直观的潜在脱靶基因/转录本列表，包含风险分和主要风险原因。
- `Validation Direction`：给出后续验证建议，例如优先做 BLAST、qPCR 或片段比对。
- 排名会综合设计得分和脱靶风险，不再只按共识评分排序。

## 非冗余推荐
- 软件默认折叠高度重叠的滑窗候选，避免前几十条只相差 1 bp。
- `Cluster Size` 表示该推荐代表了多少个相似候选。
- `Raw` 是原始滑窗候选数，`Recommended` 是折叠后的非冗余推荐数。
- 导出的 CSV/Excel/FASTA 会包含 cluster 信息，便于追踪替代窗口。

## 实验验证报告和引物
- `Export Report` 生成 Excel 报告，包含推荐理由、潜在脱靶基因/转录本、风险原因、匹配片段、T7 引物和方法说明。
- `Export Primers` 生成 T7 引物、sgRNA 克隆寡核苷酸和分型 PCR 引物订购 CSV。
- Long dsRNA 结果会自动生成普通 PCR primer 和 T7 promoter primer；订购前仍建议复核扩增唯一性。
- sgRNA 模式建议提交 CDS 序列；如果选择 mRNA/cDNA，软件会尝试推断最长 ATG-to-stop CDS，并优先推荐前段 CDS 候选。
- sgRNA 结果会显示 SpCas9 20nt spacer、NGG PAM、strand、cut site、CDS 位置、原始输入坐标、Top off-target、BbsI 兼容克隆 oligo 和切点附近 genotyping PCR primer；高风险位点建议做靶向扩增测序/ICE/TIDE/Sanger 复核。
- sgRNA 脱靶评分会先把参考/背景 FASTA 的 NRG PAM 位点索引一次，再批量评分所有候选；同一轮设计不会对 278 个候选重复全库扫描 278 次，同一软件会话里换不同基因可复用同一转录组索引。
- sgRNA off-target 只覆盖当前加载的参考/背景 FASTA；如果只加载转录组，不覆盖 intron/intergenic 区域。做 Cas9 基因组级脱靶筛查时，请加载 genome FASTA 作为参考或额外背景。
- `Save Project` / `Open Project` 可保存或恢复 `.dsforge_project` 项目文件。

## 常见问题
- 双击没有窗口：把发布包解压到英文或中文都可以的普通用户目录，例如桌面；不要直接在压缩包里运行。
- 导出 Excel 失败：请优先用发布包内的 exe，不要用源码目录里的 Python 脚本运行。
- 结果为 0 个候选：软件会弹出可能原因和建议；通常先检查目标序列长度，再切换到 `Relaxed` 模式或放宽高级设置里的 GC 范围。
