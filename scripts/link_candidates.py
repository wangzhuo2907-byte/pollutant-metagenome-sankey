"""将 PUA 候选科、属、种按基因 ID 连接到 CAZy level2 与 EC。

仅依赖 Python 标准库。.xls 输入实际为 UTF-8 制表符文本。
候选来自相关性表，不重新检验；按完整分类谱系匹配，不做名称模糊匹配。
"""
import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

SAMPLES = []
PREFIXES = ['k__','p__','c__','o__','f__','g__','s__']
RANK_END = {'Family': 5, 'Genus': 6, 'Species': 7}


def read_rows(path, required, delimiter='\t'):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        rows = csv.DictReader(handle, delimiter=delimiter)
        absent = set(required) - set(rows.fieldnames or [])
        if absent:
            raise ValueError(f'{path} 缺少字段：{sorted(absent)}')
        yield from rows


def parse_taxonomy(text):
    parts = {part.strip()[:3]: part.strip()[3:].strip()
             for part in text.split(';') if '__' in part}
    if 'k__' not in parts and 'd__' in parts:
        parts['k__'] = parts['d__']
    return parts


def extract_ec(description):
    numbers = set()
    for block in re.findall(r'\(\s*EC\s+([^)]*)\)', description, re.I):
        numbers.update(re.findall(
            r'(?<![\d.])[1-7]\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)(?![\d.])', block))
    return tuple(sorted(numbers))


def make_writer(handle, columns):
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
    writer.writeheader()
    return writer


def dump_tsv(path, columns, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = make_writer(handle, columns)
        writer.writerows(rows)


def main():
    global SAMPLES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--correlations', type=Path, required=True)
    parser.add_argument('--bacteria', type=Path, required=True)
    parser.add_argument('--annotation', type=Path, required=True)
    parser.add_argument('--level2', type=Path, required=True)
    parser.add_argument('--ec', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threshold', type=float, default=0.1, help='平均相对丰度阈值，百分数；0.1 表示 0.1%%')
    parser.add_argument('--rho', type=float, default=0.4)
    parser.add_argument('--p', type=float, default=0.05)
    parser.add_argument('--positive-only', action='store_true', help='只保留正相关；默认同时保留正负相关')
    parser.add_argument('--samples', required=True, help='样本列名，以逗号分隔')
    parser.add_argument('--pollutant-column', default=None, help='相关表的污染物列，默认识别Pollutant或PUA')
    args = parser.parse_args()
    SAMPLES = [x.strip() for x in args.samples.split(',') if x.strip()]
    if not SAMPLES or len(SAMPLES) != len(set(SAMPLES)):
        raise ValueError('样本列为空或重复')
    if not (args.threshold >= 0 and 0 <= args.rho <= 1 and 0 < args.p <= 1):
        raise ValueError('筛选阈值无效')
    for name in ['correlations','bacteria','annotation','level2','ec']:
        if not getattr(args, name).is_file():
            raise FileNotFoundError(getattr(args, name))
    args.output.mkdir(parents=True, exist_ok=True)

    # 每个候选用 (分类层级, 完整谱系) 作为唯一键，避免同名或不同层级混淆。
    candidates, significant_rows, seen_pairs = {}, [], set()
    for row in read_rows(args.correlations,
        ['Level','Taxon_key','Taxon','Rho','P_value','Status','Threshold_percent'], ','):
        pollutant_column = args.pollutant_column or ('Pollutant' if 'Pollutant' in row else 'PUA')
        if pollutant_column not in row or not row[pollutant_column].strip():
            raise ValueError('缺少污染物名称列：' + pollutant_column)
        row = dict(row)
        row['PUA'] = row[pollutant_column].strip()
        threshold = float(row['Threshold_percent'])
        if not math.isfinite(threshold) or threshold < 0:
            raise ValueError('相关表的平均丰度阈值无效')
        if row['Level'] not in RANK_END or abs(threshold-args.threshold)>1e-9:
            continue
        pair = (row['Level'], row['Taxon_key'], row['PUA'])
        if pair in seen_pairs:
            raise ValueError(f'相关表出现重复组合：{pair}')
        seen_pairs.add(pair)
        if row['Status'] != 'OK' or not row['Rho'] or not row['P_value']:
            continue
        rho, p = float(row['Rho']), float(row['P_value'])
        strength = rho if args.positive_only else abs(rho)
        if not (-1 <= rho <= 1 and 0 <= p <= 1):
            raise ValueError('相关系数或p值无效')
        if strength < args.rho or p >= args.p:
            continue
        key = (row['Level'], row['Taxon_key'])
        candidate = candidates.setdefault(key, {'Level': key[0], 'Taxon_key': key[1],
            'Taxon': row['Taxon'], 'PUA_links': [],
            'Mean_abundance_percent': row.get('Mean_abundance_percent','')})
        candidate['PUA_links'].append({'PUA': row['PUA'], 'Rho': rho, 'P_value': p})
        significant_rows.append(row)
    if not candidates:
        raise ValueError('没有满足条件的候选分类群')
    print('1/4 候选分类群：'+json.dumps(dict(Counter(k[0] for k in candidates)),ensure_ascii=False),flush=True)
    dump_tsv(args.output/'selected_PUA_correlations.tsv', list(significant_rows[0]), significant_rows)
    level2_catalog = {r['CAZy_Family'].strip() for r in read_rows(args.level2,['CAZy_Family'])}
    ec_catalog = {r['EC_ID'].strip() for r in read_rows(args.ec,['EC_ID'])}

    memberships = {}
    totals, annotated_ids = defaultdict(set), defaultdict(set)
    counts = Counter()
    # 第一遍细菌表：只保存候选基因的 ID 与分类群归属，不保存数百万行丰度。
    print('2/4 从细菌丰度表提取候选分类群的基因 ID……',flush=True)
    for row in read_rows(args.bacteria,['Reference_ID','Taxonomy']+SAMPLES):
        counts['bacterial_input_rows'] += 1
        for sample in SAMPLES:
            try:
                value = float(row[sample])
            except (ValueError, TypeError) as exc:
                raise ValueError(f'基因 {row["Reference_ID"]} 的样本 {sample} 丰度不是有效数值') from exc
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'基因 {row["Reference_ID"]} 的样本 {sample} 丰度缺失、非有限或为负')
        parts = parse_taxonomy(row['Taxonomy'])
        if parts.get('k__') != 'Bacteria':
            raise ValueError('细菌输入存在非 k__Bacteria 记录')
        member_keys = []
        for level, end in RANK_END.items():
            # 与前期 R 脚本 paste() 的 NA 占位符保持一致。
            lineage = ' | '.join(parts.get(p) or 'NA' for p in PREFIXES[:end])
            key = (level,lineage)
            if key in candidates:
                member_keys.append(key)
        if member_keys:
            gene = row['Reference_ID'].strip()
            if not gene or gene in memberships:
                raise ValueError(f'候选基因 ID 缺失或重复：{gene}')
            memberships[gene] = tuple(member_keys)
            for key in member_keys:
                totals[key].add(gene)
    counts['candidate_unique_gene_ids'] = len(memberships)
    print('候选分类群涉及去重基因数：'+str(len(memberships)),flush=True)

    # 直接流式读取原始注释，不依赖此前可能未完成的 SQLite 缓存。
    # 同一基因的重复注释只保留一次；保留一基因多家族和多描述关系。
    annotations = defaultdict(set)
    print('3/4 按 Gene_ID = Reference_ID 连接 CAZy 注释……',flush=True)
    for row in read_rows(args.annotation,['Gene_ID','CAZy_Family','Description','Taxonomy']):
        counts['annotation_input_rows'] += 1
        gene = row['Gene_ID'].strip()
        if gene in memberships:
            family = row['CAZy_Family'].strip()
            if not family:
                counts['candidate_annotation_rows_without_family'] += 1
                continue
            annotations[gene].add((family,row['Description'].strip(),row['Taxonomy'].strip()))
        if counts['annotation_input_rows'] % 500000 == 0:
            print('已读取注释行数：'+str(counts['annotation_input_rows']),flush=True)
    counts['matched_candidate_unique_gene_ids'] = len(annotations)

    families, ec_labels, ec_number_sets = defaultdict(set), defaultdict(set), defaultdict(set)
    link_columns = ['Level','Taxon_key','Taxon','Reference_ID','CAZy_level2','In_level2_catalog',
        'EC_ID','EC_numbers','In_ec_catalog','Description','Annotation_Taxonomy','PUA_links_JSON']
    print('4/4 输出科、属、种与 level2、EC 的对应结果……',flush=True)
    with open(args.output/'candidate_gene_CAZy_links.tsv','w',encoding='utf-8-sig',newline='') as handle:
        writer = make_writer(handle,link_columns)
        for gene in sorted(annotations):
            for key in memberships[gene]:
                annotated_ids[key].add(gene)
                for family,desc,tax in sorted(annotations[gene]):
                    numbers = extract_ec(desc)
                    families[(key,family)].add(gene)
                    if desc in ec_catalog or numbers:
                        ec_labels[(key,desc)].add(gene)
                    for number in numbers:
                        ec_number_sets[(key,number)].add(gene)
                    writer.writerow({'Level':key[0],'Taxon_key':key[1],'Taxon':candidates[key]['Taxon'],
                        'Reference_ID':gene,'CAZy_level2':family,'In_level2_catalog':family in level2_catalog,
                        'EC_ID':desc if desc in ec_catalog else '', 'EC_numbers':';'.join(numbers),
                        'In_ec_catalog':desc in ec_catalog,'Description':desc,'Annotation_Taxonomy':tax,
                        'PUA_links_JSON':json.dumps(candidates[key]['PUA_links'],ensure_ascii=False)})
                    counts['candidate_gene_annotation_links_across_levels'] += 1

    # 每层级分别汇总，计数单位是去重基因 ID，而不是注释行数。
    for level in RANK_END:
        summary_rows, family_rows, ec_rows, number_rows = [], [], [], []
        for key in sorted(candidates):
            if key[0] != level:
                continue
            info = candidates[key]
            base = {k:info[k] for k in ['Level','Taxon_key','Taxon']}
            f = sorted(family for (k,family) in families if k==key)
            e = sorted(desc for (k,desc) in ec_labels if k==key)
            summary_rows.append(dict(base,Mean_abundance_percent=info['Mean_abundance_percent'],
                PUA_links_JSON=json.dumps(info['PUA_links'],ensure_ascii=False),
                Total_gene_count=len(totals[key]),CAZy_gene_count=len(annotated_ids[key]),
                CAZy_gene_fraction=len(annotated_ids[key])/len(totals[key]) if totals[key] else '',
                CAZy_family_count=len(f),CAZy_families=';'.join(f),EC_label_count=len(e)))
            for family in f:
                family_rows.append(dict(base,CAZy_level2=family,In_level2_catalog=family in level2_catalog,
                    Unique_gene_count=len(families[(key,family)])))
            for desc in e:
                ec_rows.append(dict(base,EC_ID=desc if desc in ec_catalog else '',Description=desc,
                    EC_numbers=';'.join(extract_ec(desc)),In_ec_catalog=desc in ec_catalog,
                    Unique_gene_count=len(ec_labels[(key,desc)])))
            for (k,number),ids in sorted(ec_number_sets.items()):
                if k==key:
                    number_rows.append(dict(base,EC_number=number,Unique_gene_count=len(ids)))
        common = ['Level','Taxon_key','Taxon']
        dump_tsv(args.output/f'{level}_summary.tsv',common+['Mean_abundance_percent','PUA_links_JSON',
            'Total_gene_count','CAZy_gene_count','CAZy_gene_fraction','CAZy_family_count','CAZy_families','EC_label_count'],summary_rows)
        dump_tsv(args.output/f'{level}_level2.tsv',common+['CAZy_level2','In_level2_catalog','Unique_gene_count'],family_rows)
        dump_tsv(args.output/f'{level}_ec_labels.tsv',common+['EC_ID','Description','EC_numbers','In_ec_catalog','Unique_gene_count'],ec_rows)
        dump_tsv(args.output/f'{level}_ec_numbers.tsv',common+['EC_number','Unique_gene_count'],number_rows)

    # 第二遍细菌表：只输出获得 CAZy 注释的候选基因。丰度每个基因仅出现一次。
    gene_columns = ['Reference_ID']+SAMPLES+['Taxonomy','Candidate_memberships_JSON','CAZy_families',
        'EC_numbers','EC_ID_labels_JSON','Annotation_taxonomy_agrees']
    with open(args.output/'candidate_CAZy_gene_abundance.tsv','w',encoding='utf-8-sig',newline='') as handle:
        writer=make_writer(handle,gene_columns)
        for row in read_rows(args.bacteria,['Reference_ID','Taxonomy']+SAMPLES):
            gene=row['Reference_ID'].strip()
            if gene not in annotations:
                continue
            links=annotations[gene]
            agreement=all(parse_taxonomy(tax)==parse_taxonomy(row['Taxonomy']) for _,_,tax in links)
            counts['matched_genes_with_taxonomy_disagreement'] += not agreement
            writer.writerow(dict({k:row[k] for k in ['Reference_ID']+SAMPLES+['Taxonomy']},
                Candidate_memberships_JSON=json.dumps([{'Level':k[0],'Taxon_key':k[1]} for k in memberships[gene]],ensure_ascii=False),
                CAZy_families=';'.join(sorted({f for f,_,_ in links})),
                EC_numbers=';'.join(sorted({n for _,d,_ in links for n in extract_ec(d)})),
                EC_ID_labels_JSON=json.dumps(sorted({d for _,d,_ in links if d in ec_catalog}),ensure_ascii=False),
                Annotation_taxonomy_agrees=agreement))
            counts['output_abundance_gene_rows'] += 1
    if counts['output_abundance_gene_rows'] != len(annotations):
        raise ValueError('输出基因数量核对失败')
    report=dict(counts)
    report['genes_with_multiple_EC_labels']=sum(
        len({d for _,d,_ in links if d in ec_catalog})>1 for links in annotations.values())
    report['EC_interpretation']='EC 关系继承源文件 Description；同一基因/家族的多种酶描述需核查上游注释规则，不能视为已确认的全部酶活。'
    report['candidate_counts_by_level']=dict(Counter(k[0] for k in candidates))
    report['candidate_CAZy_gene_counts_by_level']={level:sum(len(annotated_ids[k]) for k in candidates if k[0]==level) for level in RANK_END}
    report['families_absent_from_level2_catalog']=sorted({f for (_,f) in families if f not in level2_catalog})
    report['ec_descriptions_absent_from_ec_catalog']=sorted({d for (_,d) in ec_labels if d not in ec_catalog})
    report['samples']=SAMPLES
    report['settings']={'threshold_percent':args.threshold,'rho':args.rho,'p':args.p,'positive_only':args.positive_only}
    report['inputs']={k:str(getattr(args,k)) for k in ['correlations','bacteria','annotation','level2','ec']}
    (args.output/'matching_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    (args.output/'README.txt').write_text(
        '候选按平均相对丰度阈值和至少一种 PUA 的相关条件筛选，正负相关默认均保留。\n'
        '基因归属使用细菌表的完整分类谱系；CAZy 使用 Gene_ID=Reference_ID 连接。\n'
        'Family/Genus/Species_summary：候选的 CAZy 覆盖、家族清单和 PUA 关联。\n'
        '*_level2：分类群与 CAZy 家族对应及去重基因数；*_ec_labels：EC 原始描述类别；*_ec_numbers：EC 编号。\n'
        'EC_ID 仅在完整描述与 EC 汇总表精确匹配时填写；未匹配描述仍保留提取出的 EC 编号。\n'
        'EC 关系继承源文件 Description；同一基因/家族的多种酶描述需核查上游注释规则，不能视为已确认的全部酶活。\n'
        'candidate_gene_CAZy_links 是关系明细，不携带丰度；candidate_CAZy_gene_abundance 每个基因一行。\n'
        'CAZy_gene_fraction 是基因记录数比例，不是丰度比例。\n'
        '同一基因可有多个家族/EC，也可属于候选科和属；各类群/功能基因数不可直接相加作为去重总数。\n'
        'level2/EC 汇总表仅用于类别核对，不将家族总丰度复制到基因。\n'
        '未匹配到 CAZy 不等于该分类群完全没有该功能；结果不代表表达、酶活或 PUA 的因果作用。\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
    print('完成：'+str(args.output),flush=True)


if __name__=='__main__':
    main()
