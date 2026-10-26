# 输入与命令

## 规范化输入

所有基因ID保持字符串，不截断、不去掉样本前缀。TSV/CSV用UTF-8；真实Excel可在上游读取或先导出。每种数据均明确行粒度。

| 表 | 必需列 | 规则 |
|---|---|---|
| 全细菌基因丰度TSV | `Reference_ID, Taxonomy, 样本列…` | 每ID一行，样本列数可变。注释基因丰度不是细胞绝对数量。不要只输入CAZy子集计算全细菌相对丰度。 |
| 污染物浓度Excel/CSV/TSV | `SampleID, 指定污染物列…` | 上游R可识别Station前缀及中文说明；样本ID需与丰度列精确一致。其他情况先规范化。Layer可选；未提供时标为Unspecified，不按字母猜测水层。 |
| 相关结果CSV | `Level, Taxon_key, Taxon, PUA或Pollutant, Rho, P_value, Status, Threshold_percent` | Level为Family/Genus/Species；阈值用百分数0.1表示0.1%。可带N和BH q。Taxon_key为完整谱系用` | `连接，缺失层级为NA。 |
| 逐基因CAZy TSV | `Gene_ID, CAZy_Family, Description, Taxonomy` | 一基因可多行。Description中的EC需注明来源可信度。 |
| level2类别TSV | `CAZy_Family` | 其他列可含汇总丰度；仅核对类别，不用于反推基因。 |
| EC类别TSV | `EC_ID` | 本脚本处理EC_ID为完整描述的格式；若是纯编号，应先调整脚本核对逻辑，不能直接按描述连接后宣称成功。 |

`Taxon_key`示例：

```text
Family:  Bacteria | Bacteroidota | Flavobacteriia | Flavobacteriales | Crocinitomicaceae
Genus:   Bacteria | … | FamilyName | GenusName
Species: Bacteria | … | FamilyName | GenusName | SpeciesName
```

若已有可信的候选名单但无相关p/rho，不编造这些值。联合注释可按名单提取ID；绘图第一段标为用户提供的关系，并相应改写脚本，不冒充显著关联。

## 执行

在下列命令中，`SKILL`代表本技能安装目录，所有路径由当前任务替换。

### 可选：从原始数据计算相关性

```text
Rscript SKILL/scripts/correlate_taxa.R --bacteria bacteria.tsv --pollutants concentrations.xlsx --output correlations --pollutant-columns pC7,pC8,hpC7 --permutations 9999 --min-n 8
```

脚本不假设固定样本数量；按浓度表样本名匹配丰度列。默认按0.1%入口计算，另报告0.5%和1%子集的校正结果；原始p不因丰度阈值子集改变。数据少于`min-n`不会强行降低样本门槛。零值与缺失值不同，检测下限值须先制定处理规则。

### 候选与CAZy连接

```text
python SKILL/scripts/link_candidates.py --correlations correlations/correlation_by_threshold.csv --bacteria bacteria.tsv --annotation Unigenes.CAZy.tax.xls --level2 Unigenes.absolute.level2.xls --ec Unigenes.absolute.ec.xls --output matched --threshold 0.1 --rho 0.4 --p 0.05 --samples B12,D12,B10,D10,B3,B4,B2,B1,B9,D9,B11
```

`--samples`明确列出本次样本；不省略也不默认使用案例样本。默认接受相关表的`PUA`或`Pollutant`列，可用`--pollutant-column`指定其他列名。`--positive-only`只筛正相关。输出内部字段`PUA_links_JSON`可装任意污染物名称。

### 桑基图与交互版

```text
python SKILL/scripts/sankey.py --input matched --output figures --top-families 12 --top-ec 15 --pollutant-label Pollutants
python SKILL/scripts/interactive.py --input figures/sankey_display_data.json --output figures/interactive_sankey.html
```

PNG供预览，SVG/PDF供编辑排版。完整四元路径在`full_annotation_paths.tsv`；显示类别映射在`display_category_mapping.tsv`。README与审计JSON保留算法口径、参数和缺失类别。

## 本次PUA案例（仅作为示例）

六项指标pC7/pC8/pC10/hpC7/hpC8/hpC10；11个独立样本，B表层、D底层。候选按平均相对丰度≥0.1%、|rho|≥0.4、原始p<0.05选择。结果为5科、6属、11种级标签，并非22株独立细菌。该案例未控制水层且无校正后显著结果，只作探索。新任务不得直接套用这些样本、候选数量或结论。
