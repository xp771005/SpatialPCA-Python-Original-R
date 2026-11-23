# SpatialPCA Python 与原版 R 对照

本次已实际运行作者未修改的 R 模型核心，与独立 Python 实现在同一个 Y 和 K 上进行 9 组模拟对照，全部通过数值门槛。历史本地结果与日志见 results/local_reference；GitHub Actions #2 也已实际通过 18 项测试和 9 组对照，运行记录见仓库根目录 experiment_results.md。

## 需要更正的旧结论

旧版本把“去掉 log(S)”的 Python 目标函数误标为原版 R 表达式。作者源码固定版本 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66 明确包含 log(S)。旧结论撤回。对应 Python 模式改名为 no_log_ablation，只用于人为改变目标函数的消融实验。

R Runner 的括号语法错误及 Z 矩阵错误转置已修复。原版 R 的 SpatialPCs 已经是 d×n，不应再次转置。原作者源码未修改。

## 与 NJIT 任务的关系

Zhi Wei 老师要求理解模型、补全跳步、核对笔误和实现 Python。现在核心实现及模拟跨语言数值对照已有证据。推导笔记保留了对 B 的积分、Z 积分、噪声方差剖面估计、tau 优化和核低秩处理的中间步骤。

本轮并未重新逐页核对论文和所有附录，所以不能把先前记录的正文维度或附录笔误说成这次重新确认的结论。特别需要把论文、博士论文附录和某个固定版本 R 代码区分开。

## 当前状态与范围

仓库已公开，源码可直接浏览。当前 Actions 直接运行根目录源码，已移除历史 ZIP 解压和私有仓库检查。LaTeX 报告与 PDF 分别位于 docs/SpatialPCA_Short_Note.tex 和 docs/SpatialPCA_Short_Note.pdf。当前结果仅覆盖模拟数据、截距项和完整核；不覆盖 SCTransform/SPARK、真实组织数据和全部论文图表。
