"""Build three conserved PUA–taxon–CAZy–EC Sankey diagrams.

Requires matplotlib for PNG/SVG/PDF; all joins/aggregation use the standard library.
Fractional equal allocation is a display convention, not measured abundance or causality.
"""
import argparse
import csv
import json
import math
from collections import defaultdict, Counter
from pathlib import Path

PUAS = []
LEVELS = ['Species','Family','Genus']
OTHER_FAMILY = 'Other CAZy families'
OTHER_EC = 'Other EC numbers'
NO_EC = 'EC unassigned'
COLORS = {'GH':'#718CAC','GT':'#80A594','CBM':'#A190B1','CE':'#BD9E72','PL':'#7AA4AC','AA':'#AF8A89','Other':'#A3A8B0'}


def read_tsv(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        yield from csv.DictReader(f,delimiter='\t')


def write_tsv(path,fields,rows):
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t'); w.writeheader(); w.writerows(rows)


def family_color(label):
    return next((color for prefix,color in COLORS.items() if label.startswith(prefix)),COLORS['Other'])


def draw_static(level,paths,taxa,total,out,settings,top_families,top_ec,pollutant_label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch,Rectangle
    plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
    edges=defaultdict(float); edge_meta={}
    weights=defaultdict(float); labels={}; layer={}
    for r in paths:
        nodes=[('PUA',r['pua']),('Taxon',r['taxon_key']),('Family',r['family']),('EC',r['ec'])]
        for i,node in enumerate(nodes):
            labels[node]=[r['pua'],r['taxon'],r['family'],r['ec']][i]
            layer[node]=i
        for a,b in zip(nodes,nodes[1:]):
            edges[(a,b)]+=r['weight']
            if a[0]=='PUA': edge_meta[(a,b)]=(r['rho'],r['p'])
    for (a,b),v in edges.items():
        if layer[a]==0: weights[a]+=v
        weights[b]+=v
    # Show all configured pollutants, including sources with no association at this rank.
    for p in PUAS:
        node=('PUA',p);labels[node]=p;layer[node]=0;weights.setdefault(node,0.)
    for node in weights:
        if layer[node] in (1,2):
            outgoing=sum(v for (a,b),v in edges.items() if a==node)
            if not math.isclose(outgoing,weights[node],rel_tol=1e-9,abs_tol=1e-7):
                raise AssertionError('Nonconserved node '+str(node))
    columns=[]
    for col in range(4):
        nodes=[n for n in weights if layer[n]==col]
        nodes.sort(key=(lambda n:PUAS.index(n[1])) if col==0 else (lambda n:(-weights[n],labels[n])))
        columns.append(nodes)
    W,H=2500,1420
    if total <= 0:
        fig,ax=plt.subplots(figsize=(12,6));ax.axis('off')
        ax.text(.5,.55,level+' — no matched candidate CAZy genes',ha='center',va='center')
        for ext in ['png','svg','pdf']:fig.savefig(out/('sankey_'+level.lower()+'.'+ext),dpi=180)
        plt.close(fig);return
    y0,gap=235,23
    max_n=max(map(len,columns))
    chart_h=max(1010,gap*(max_n-1)+600)
    H=chart_h+410
    scale=(chart_h-gap*(max_n-1)-6)/total
    xs=[185,850,1500,2070]; width=17
    pos={}
    for col,nodes in enumerate(columns):
        occupied=sum(max(weights[n]*scale,1.5) for n in nodes)+gap*(len(nodes)-1)
        y=y0+(chart_h-occupied)/2
        for n in nodes:
            height=max(weights[n]*scale,1.5)
            pos[n]=(xs[col],y,height)
            y+=height+gap
    starts={};ends={}
    for node in weights:
        outgoing=sorted(((b,v) for (a,b),v in edges.items() if a==node),key=lambda z:pos[z[0]][1])
        incoming=sorted(((a,v) for (a,b),v in edges.items() if b==node),key=lambda z:pos[z[0]][1])
        offset=0
        for b,v in outgoing:starts[(node,b)]=pos[node][1]+offset;offset+=v*scale
        offset=0
        for a,v in incoming:ends[(a,node)]=pos[node][1]+offset;offset+=v*scale
    fig=plt.figure(figsize=(22,12.5*H/1420),facecolor='white')
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,W);ax.set_ylim(H,0);ax.axis('off')
    title=f'{pollutant_label} – candidate {level.lower()} – CAZy family – EC'
    ax.text(70,62,title,fontsize=22,weight='bold',color='#223340')
    ax.text(70,111,f"{len(taxa)} matched candidate taxa  |  {total:,} unique CAZy genes  |  mean abundance ≥{settings['threshold_percent']}%, |rho| ≥{settings['rho']}, raw p <{settings['p']}",fontsize=12,color='#5B6570')
    ax.plot([75,110],[151,151],lw=6,color='#C77057',solid_capstyle='round')
    ax.text(122,157,'Positive association',fontsize=11,color='#3D4851')
    ax.plot([540,575],[151,151],lw=6,color='#408B99',solid_capstyle='round')
    ax.text(587,157,'Negative association',fontsize=11,color='#3D4851')
    ax.text(1090,157,'Downstream ribbons: inherited gene annotations',fontsize=11,color='#5B6570')
    headings=[pollutant_label,{'Family':'Bacterial family','Genus':'Bacterial genus','Species':'Species-level label'}[level],'CAZy family (level2)','EC number']
    for x,heading in zip(xs,headings):
        ax.text(x,205,heading,fontsize=13,weight='bold',color='#263A49',ha='center')
    for (a,b),value in sorted(edges.items(),key=lambda kv:-kv[1]):
        x1=pos[a][0]+width;x2=pos[b][0]
        y1=starts[(a,b)];y2=ends[(a,b)];h=value*scale
        bend=(x2-x1)*.48
        points=[(x1,y1),(x1+bend,y1),(x2-bend,y2),(x2,y2),(x2,y2+h),
                (x2-bend,y2+h),(x1+bend,y1+h),(x1,y1+h),(x1,y1)]
        codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,
               MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
        if layer[a]==0: color='#C77057' if edge_meta[(a,b)][0]>0 else '#408B99';alpha=.38
        elif layer[a]==1: color=family_color(labels[b]);alpha=.30
        else: color=family_color(labels[a]);alpha=.26
        ax.add_patch(PathPatch(MPath(points,codes),facecolor=color,edgecolor='none',alpha=alpha))
    for node,(x,y,h) in pos.items():
        col=layer[node]
        color=['#415B70','#546D78',family_color(labels[node]),'#77818B'][col]
        if weights[node]==0:color='#D6DADF'
        ax.add_patch(Rectangle((x,y),width,h,facecolor=color,edgecolor='none'))
        if col in (0,1): tx=x-12;align='right'
        else: tx=x+width+12;align='left'
        ax.text(tx,y+h/2,labels[node],va='center',ha=align,fontsize=10.7 if col==1 else 11.2,
            color='#23323C',bbox=dict(facecolor='white',alpha=.84,pad=1.4,edgecolor='none'))
    ax.text(70,H-116,'Ribbon width = fractionally allocated gene annotations; each gene contributes a total weight of 1 within each diagram.',fontsize=11,color='#4C5962')
    ax.text(70,H-80,f'Top {top_families} CAZy families and top {top_ec} EC numbers shown; other categories retained in Other. Unassigned EC retained.',fontsize=10.5,color='#5C6670')
    ax.text(70,H-45,'Exploratory, not causal; EC are inherited candidates. Unmatched taxa remain in input summaries; merged nodes do not prove every full path.',fontsize=10.5,color='#5C6670')
    stem=out/f'sankey_{level.lower()}'
    for ext in ['png','svg','pdf']:
        fig.savefig(str(stem)+'.'+ext,dpi=180,facecolor='white')
    plt.close(fig)


def main():
    global PUAS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--top-families',type=int,default=12)
    parser.add_argument('--top-ec',type=int,default=15)
    parser.add_argument('--pollutant-label',default='Pollutants')
    parser.add_argument('--pollutant-order',default=None)
    parser.add_argument('--threshold',type=float,default=.1)
    parser.add_argument('--rho',type=float,default=.4)
    parser.add_argument('--p',type=float,default=.05)
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    if args.top_families < 1 or args.top_ec < 1:raise ValueError('Top counts must be positive')
    settings={'threshold_percent':args.threshold,'rho':args.rho,'p':args.p}
    report_path=args.input/'matching_report.json'
    if report_path.exists():settings.update(json.loads(report_path.read_text(encoding='utf-8')).get('settings',{}))
    associations={};taxa_info={};observed_pollutants=[]
    for r in read_tsv(args.input/'selected_PUA_correlations.tsv'):
        key=(r['Level'],r['Taxon_key'])
        taxa_info[key]=r['Taxon']
        if r['PUA'] not in observed_pollutants:observed_pollutants.append(r['PUA'])
        associations.setdefault(key,{})[r['PUA']]=(float(r['Rho']),float(r['P_value']))
    PUAS = [x.strip() for x in args.pollutant_order.split(',')] if args.pollutant_order else observed_pollutants
    if len(PUAS)!=len(set(PUAS)) or not set(observed_pollutants).issubset(PUAS):raise ValueError('污染物顺序缺失或重复')
    # Collapse repeated enzyme descriptions before assigning gene weight.
    gene_families=defaultdict(lambda:defaultdict(set))
    family_present=set();no_ec_present=set()
    for r in read_tsv(args.input/'candidate_gene_CAZy_links.tsv'):
        key=(r['Level'],r['Taxon_key'],r['Reference_ID'])
        f=r['CAZy_level2'];family_present.add((key,f))
        if r['EC_numbers'].strip():
            gene_families[key][f].update(x.strip() for x in r['EC_numbers'].split(';') if x.strip())
        else:
            no_ec_present.add((key,f));gene_families[key][f]
    full_paths=[];payload={'meta':{'unit':'fractionally allocated gene annotations','top_families':args.top_families,'top_ec':args.top_ec,
        'pollutant_order':PUAS,'pollutant_label':args.pollutant_label,'threshold_percent':settings['threshold_percent'],'rho_cutoff':settings['rho'],'p_cutoff':settings['p']},'levels':{}}
    summaries={};grouping_rows=[]
    for level in LEVELS:
        keys_at_level=[key for key in taxa_info if key[0]==level]
        display_names={}
        for _,taxon_key in keys_at_level:
            raw=taxa_info[(level,taxon_key)]
            peers=[k for lev,k in keys_at_level if taxa_info[(lev,k)]==raw]
            label=raw
            if len(peers)>1:
                parts=taxon_key.split(' | ')
                for depth in range(2,len(parts)+1):
                    suffix=' / '.join(parts[-depth:])
                    if all(other==taxon_key or ' / '.join(other.split(' | ')[-depth:])!=suffix for other in peers):
                        label=suffix
                        break
                else:
                    label=taxon_key
            display_names[taxon_key]=label
        allocation=defaultdict(float);gene_ids=set();taxon_gene_ids=defaultdict(set)
        for (lev,taxon,gene),family_map in gene_families.items():
            if lev!=level:continue
            gene_ids.add(gene);taxon_gene_ids[taxon].add(gene)
            n_family=len(family_map)
            for f,ecs in family_map.items():
                # An unnumbered description is not an additional unknown enzyme when
                # the same gene-family also has numbered descriptions.
                targets=ecs or {NO_EC}
                w=1.0/(n_family*len(targets))
                for ec in targets:allocation[(taxon,f,ec)]+=w
        total=len(gene_ids)
        if not total:
            payload['levels'][level]={'taxa':[],'paths':[],'total_genes':0}
            summaries[level]={'taxa':0,'unique_genes':0,'conserved_total_weight':0,
                'candidate_taxa_before_matching':sum(lev==level for lev,t in associations)}
            draw_static(level,[],[],0,out,settings,args.top_families,args.top_ec,args.pollutant_label)
            continue
        if not math.isclose(sum(allocation.values()),total,rel_tol=1e-9):raise AssertionError('Gene allocation failed')
        family_weights=Counter();ec_weights=Counter()
        for (_,f,e),w in allocation.items():family_weights[f]+=w;ec_weights[e]+=w
        top_f={x for x,_ in sorted(family_weights.items(),key=lambda kv:(-kv[1],kv[0]))[:args.top_families]}
        top_e=set([x for x,_ in sorted(ec_weights.items(),key=lambda kv:(-kv[1],kv[0])) if x!=NO_EC][:args.top_ec])
        displayed=defaultdict(float)
        for (taxon,f,e),w in sorted(allocation.items()):
            assoc=associations[(level,taxon)]
            for pua,(rho,p) in assoc.items():
                weight=w/len(assoc)
                full_paths.append({'Level':level,'PUA':pua,'Taxon_key':taxon,'Taxon':taxa_info[(level,taxon)],
                    'CAZy_family':f,'EC_number':e,'Weight':weight,'Rho':rho,'P_value':p})
                shown_f=f if f in top_f else OTHER_FAMILY
                shown_e=e if e in top_e or e==NO_EC else OTHER_EC
                displayed[(pua,taxon,shown_f,shown_e)]+=weight
        paths=[]
        for (pua,taxon,f,e),w in sorted(displayed.items()):
            rho,p=associations[(level,taxon)][pua]
            paths.append(dict(pua=pua,taxon=display_names[taxon],taxon_key=taxon,family=f,ec=e,weight=w,rho=rho,p=p))
        if not math.isclose(sum(r['weight'] for r in paths),total,rel_tol=1e-9):raise AssertionError('Grouped allocation failed')
        taxa=[{'key':t,'label':display_names[t],'gene_count':len(ids)} for t,ids in sorted(taxon_gene_ids.items())]
        payload['levels'][level]={'taxa':taxa,'paths':paths,'total_genes':total}
        summaries[level]={'taxa':len(taxa),'candidate_taxa_before_matching':sum(lev==level for lev,t in associations),'unique_genes':total,'CAZy_families_before_grouping':len(family_weights),
            'EC_numbers_before_grouping':len(ec_weights),'display_paths':len(paths),'conserved_total_weight':sum(r['weight'] for r in paths)}
        for f,w in family_weights.items():grouping_rows.append(dict(Level=level,Stage='CAZy',Original=f,Displayed=f if f in top_f else OTHER_FAMILY,Weight=w))
        for e,w in ec_weights.items():grouping_rows.append(dict(Level=level,Stage='EC',Original=e,Displayed=e if e in top_e or e==NO_EC else OTHER_EC,Weight=w))
        print('Plotting '+level+': '+json.dumps(summaries[level]),flush=True)
        draw_static(level,paths,taxa,total,out,settings,args.top_families,args.top_ec,args.pollutant_label)
    write_tsv(out/'full_annotation_paths.tsv',['Level','PUA','Taxon_key','Taxon','CAZy_family','EC_number','Weight','Rho','P_value'],full_paths)
    write_tsv(out/'display_category_mapping.tsv',['Level','Stage','Original','Displayed','Weight'],grouping_rows)
    (out/'sankey_display_data.json').write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
    (out/'sankey_summary.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'README.txt').write_text(
        '三张桑基图：污染物 → 候选细菌科/属/种 → CAZy level2 家族 → EC 编号。\n'
        f"第一段为用户筛选的探索关联；平均丰度≥{settings['threshold_percent']}%、|rho|≥{settings['rho']}、原始p<{settings['p']}。橙红=正，青蓝=负。\n"
        '后两段来自 Reference_ID=Gene_ID 的注释连接，不是额外显著性检验或因果链。\n'
        '流量单位：分摊基因注释权重。每个基因在每张图贡献总权重1。\n'
        '具体为：在该分类群显著PUA数、该基因家族数、该基因-家族的去重EC数之间等分。\n'
        'EC缺失时保留EC unassigned；不完整编号保留横线。重复酶描述中的同一EC只计一次。\n'
        f'显示前{args.top_families}个家族和前{args.top_ec}个EC，其余合并Other但未删除；排名基于该层级分摊权重。\n'
        '这是一种制图分配约定，不是基因丰度、PUA影响强度、实际表达或酶活。\n'
        'EC是源文件继承的预测候选；多描述需核查上游注释规则。\n'
        '合并节点会失去跨层的路径身份；核对具体分类群-家族-EC请查full_annotation_paths.tsv。\n'
        '科、属、种独立作图，不能相加为独立总数；种级标签含暂定分类。\n'
        'PNG用于预览，SVG/PDF用于编辑和排版；HTML可筛选单个分类群。\n',encoding='utf-8')
    print('DONE '+str(out),flush=True)


if __name__=='__main__':main()
