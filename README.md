# 污染物-宏基因桑吉图skill

用于将污染物相关的细菌科、属、种候选与 CAZy level2、EC 按基因 ID 联合注释，生成分层桑基图及离线交互图。

## 能做什么

- 可选：从污染物浓度和全细菌基因丰度计算探索性相关。
- 按完整分类谱系提取候选基因，再以 `Reference_ID = Gene_ID` 连接功能注释。
- 保留多家族、多 EC、未匹配和分类冲突，避免重复计算丰度。
- 生成科、属、种三张 PNG/SVG/PDF，以及可筛选单个分类群的离线 HTML。
- 输出完整路径数据并检查分摊基因权重守恒。

## 使用

将仓库内容放入 Codex 的技能目录下的 `pollutant-metagenome-sankey` 文件夹。

调用：`$pollutant-metagenome-sankey`，或使用显示名称“污染物-宏基因桑吉图skill”。

输入字段及命令见 [输入与命令](references/data-and-commands.md)，统计与权重定义见 [统计与制图口径](references/methods.md)。

Python 联合注释和 HTML 生成仅依赖标准库；静态绘图需要 matplotlib。可选上游 R 分析需要 data.table、readxl、ggplot2。

## 解释边界

第一段是探索性关联，后两段是预测注释对应，不构成因果或代谢通路。线宽默认是等权分摊的基因注释权重，不是丰度、相关强度或酶活。EC 多重描述需核查上游注释规则。

本仓库不包含原始样本、宏基因组数据或个人分析结果。
