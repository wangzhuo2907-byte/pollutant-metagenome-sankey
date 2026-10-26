"""Build the offline four-column 污染物 / taxon / CAZy / EC explorer."""
from pathlib import Path
import json
import argparse


HTML = r'''<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>污染物—细菌—CAZy—EC 注释关系</title>
<style>
:root{color-scheme:light;--ink:#173047;--muted:#52677a;--line:#d8e2e9;--bg:#f4f7fa;--blue:#365f7a;--pos:#c97750;--neg:#527ba6}
*{box-sizing:border-box}body{margin:0;background:var(--bg);font-family:"Segoe UI","Microsoft YaHei",Arial,sans-serif;color:var(--ink)}
main{max-width:1600px;margin:auto;padding:28px 30px 34px}header{display:flex;gap:24px;align-items:flex-start;justify-content:space-between}
.eyebrow{font-size:12px;font-weight:600;letter-spacing:.16em;color:#477589;margin-bottom:9px}h1{font-size:28px;line-height:1.3;margin:0 0 9px;font-weight:650}p{margin:0;line-height:1.65}.intro{font-size:14px;color:var(--muted);max-width:850px}
.badge{font-size:12px;border:1px solid var(--line);border-radius:18px;padding:8px 12px;background:white;white-space:nowrap;margin-top:5px}
.toolbar{margin:23px 0 14px;display:flex;align-items:center;justify-content:space-between;gap:20px;flex-wrap:wrap}
.tabs{display:flex;background:#e6edf2;padding:4px;border-radius:8px;gap:3px}.tab{border:0;background:transparent;color:var(--muted);font:inherit;padding:9px 24px;border-radius:5px;cursor:pointer;font-size:14px}.tab[aria-selected=true]{background:white;color:var(--ink);box-shadow:0 1px 3px #1232;font-weight:600}
.filter{display:flex;gap:10px;align-items:center;font-size:14px}.filter select{font:inherit;color:var(--ink);background:white;border:1px solid #bfcfda;border-radius:6px;padding:8px 32px 8px 10px;max-width:480px}
.panel{border:1px solid var(--line);background:white;border-radius:10px;overflow:hidden}.panel-head{display:flex;justify-content:space-between;gap:18px;align-items:center;padding:17px 21px;border-bottom:1px solid #e8edf1;flex-wrap:wrap}
.metrics{font-size:13px;color:var(--muted)}.metrics strong{font-weight:650;color:var(--ink)}.legend{display:flex;gap:18px;font-size:12px;color:var(--muted);flex-wrap:wrap}.legend span{display:flex;align-items:center;gap:7px}.swatch{width:20px;height:5px;border-radius:2px;background:var(--pos)}.swatch.neg{background:var(--neg)}.swatch.neutral{background:#9bb0bd}
.scroll{overflow-x:auto;padding:0 10px}.scroll svg{display:block;width:100%;min-width:1150px;height:auto}.column-title{font-size:14px;font-weight:600;fill:#284b60}.column-sub{font-size:11px;fill:#738493}.node-label{font-size:12px;fill:#1e3b4e;pointer-events:none}.node rect{stroke:white;stroke-width:.5px;cursor:pointer}.node:focus rect{stroke:#142f44;stroke-width:2px}.node-value{font-size:10px;fill:#778b98}.link{cursor:pointer;transition:opacity .14s;stroke:none}.node{outline:none}.guide{font-size:12px;color:var(--muted);padding:11px 21px 15px;border-top:1px solid #edf1f4;min-height:45px}
.notes{margin-top:17px;display:grid;grid-template-columns:1fr 1fr;gap:14px}.note{padding:15px 18px;background:#eaf0f4;border-radius:7px;font-size:12px;color:#4b6173}.note strong{display:block;color:#294a5f;margin-bottom:5px;font-size:13px}.note p+p{margin-top:6px}.empty{font-size:16px;fill:#63788a}
@media(max-width:760px){main{padding:20px 14px}header{display:block}.badge{display:inline-block}h1{font-size:23px}.notes{grid-template-columns:1fr}.filter{width:100%;flex-wrap:wrap}.filter select{max-width:100%;width:100%}.tab{padding:8px 20px}.panel-head{padding:14px}.legend{gap:12px}}
</style>
</head>
<body>
<main>
<header><div><div class="eyebrow">METAGENOMIC ANNOTATION EXPLORER</div><h1>污染物—细菌—CAZy—EC 注释关系</h1><p class="intro">查看与污染物显著相关且匹配到CAZy基因的候选分类群，以及由基因 ID 连接得到的家族和 EC 注释。未匹配分类群保留在输入摘要中。科、属、种独立展示。</p></div><div class="badge" id="settings-badge">宏基因组 · 探索性关联</div></header>
<div class="toolbar"><div class="tabs" role="tablist" aria-label="分类层级"><button class="tab" type="button" role="tab" data-level="Family" aria-selected="true">科 Family</button><button class="tab" type="button" role="tab" data-level="Genus" aria-selected="false">属 Genus</button><button class="tab" type="button" role="tab" data-level="Species" aria-selected="false">种 Species</button></div><label class="filter">查看分类群 <select id="taxon-filter" aria-label="筛选候选分类群"></select></label></div>
<section class="panel" aria-label="桑基图"><div class="panel-head"><div class="metrics" id="metrics"></div><div class="legend"><span><i class="swatch"></i>正相关</span><span><i class="swatch neg"></i>负相关</span><span><i class="swatch neutral"></i>基因注释对应</span></div></div><div class="scroll"><svg id="chart" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="污染物、细菌分类群、CAZy 家族与 EC 的四列桑基图"></svg></div><div class="guide" id="selection" role="status">悬停查看数值；点击节点突出相邻连接，再次点击或点击空白处取消。</div></section>
<div class="notes"><div class="note"><strong>线宽如何解释</strong><p>线宽为等权分摊后的基因注释当量，不是相关系数、表达量或物质通量。每个基因的权重在该分类群关联的 污染物、该基因的 CAZy 家族及相应 EC 之间分摊，每个层级内总权重守恒。</p><p id="group-note"></p></div><div class="note"><strong>关系的证据范围</strong><p>第一段为用户设定筛选条件下的候选关联（参数见图上方）；协变量控制以输入分析为准，原始 p 不表示 FDR 校正后显著。后两段继承源文件的功能注释，多重 EC 不等于已确认的多种酶活。</p><p>流向用于展示对应关系，不代表因果或代谢通路。高亮仅显示相邻连接；合并节点不能证明任意跨四列路径都成立。</p></div></div>
</main>
<script id="sankey-data" type="application/json">__DATA__</script>
<script>
'use strict';
const DATA=JSON.parse(document.getElementById('sankey-data').textContent);
document.getElementById('settings-badge').textContent=`均值≥${DATA.meta.threshold_percent}% · |ρ|≥${DATA.meta.rho_cutoff} · 原始p<${DATA.meta.p_cutoff}`;
const svg=document.getElementById('chart'),filter=document.getElementById('taxon-filter');
const metrics=document.getElementById('metrics'),selection=document.getElementById('selection');
const NS='http://www.w3.org/2000/svg';
const POLLUTANT_ORDER=DATA.meta.pollutant_order || [...new Set(Object.values(DATA.levels).flatMap(x=>x.paths.map(p=>p.pua)))];
const LEVEL_NAMES={Family:'科',Genus:'属',Species:'种级标签'};
const nodeColors=['#527b91','#527184','#669b9a','#8292a6'];
let currentLevel='Family',selectedNode=null,view=null;
const fmt=x=>Number(x).toLocaleString('zh-CN',{maximumFractionDigits:2});
function el(tag,attrs={},parent=svg){const e=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(parent)parent.appendChild(e);return e;}
function title(parent,text){const t=el('title',{},parent);t.textContent=text;}
function text(x,y,value,attrs={},parent=svg){const t=el('text',{x,y,...attrs},parent);t.textContent=value;return t;}
function escapeHTML(x){return String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function safeP(x){return Number.isFinite(Number(x))?Number(x).toPrecision(3):'未提供';}
function key(stage,label){return stage+'|'+label;}
function prepare(paths){
 const nodes=new Map(),edgeMap=new Map();
 function node(stage,id,label){const k=key(stage,id);if(!nodes.has(k))nodes.set(k,{id:k,stage,label,in:0,out:0,linksIn:[],linksOut:[]});return nodes.get(k);}
 for(const path of paths){
  const w=Number(path.weight);if(!Number.isFinite(w)||w<=0)continue;
  const chain=[node(0,path.pua,path.pua),node(1,path.taxon_key,path.taxon),node(2,path.family,path.family),node(3,path.ec,path.ec)];
  for(let s=0;s<3;s++){
   const a=chain[s],b=chain[s+1],k=a.id+'\u001f'+b.id;
   let edge=edgeMap.get(k);if(!edge){edge={source:a,target:b,stage:s,value:0,rho:path.rho,p:path.p};edgeMap.set(k,edge);}
   edge.value+=w;
  }
 }
 const edges=[...edgeMap.values()];
 for(const e of edges){e.source.out+=e.value;e.target.in+=e.value;e.source.linksOut.push(e);e.target.linksIn.push(e);}
 const columns=[0,1,2,3].map(s=>[...nodes.values()].filter(n=>n.stage===s));
 const taxOrder=new Map(DATA.levels[currentLevel].taxa.map((t,i)=>[t.key,i]));
 columns[0].sort((a,b)=>POLLUTANT_ORDER.indexOf(a.label)-POLLUTANT_ORDER.indexOf(b.label));
 columns[1].sort((a,b)=>(taxOrder.get(a.id.slice(2))??999)-(taxOrder.get(b.id.slice(2))??999)||a.label.localeCompare(b.label));
 for(const col of columns.slice(2))col.sort((a,b)=>b.in-a.in||a.label.localeCompare(b.label));
 // Reduce crossings with weighted barycentres while keeping the 污染物 and taxa orders stable.
 for(let repeat=0;repeat<3;repeat++){
  for(let s=2;s<=3;s++){
   const previous=new Map(columns[s-1].map((n,i)=>[n.id,i]));
   const bary=n=>n.linksIn.reduce((a,e)=>a+(previous.get(e.source.id)??0)*e.value,0)/(n.in||1);
   columns[s].sort((a,b)=>bary(a)-bary(b)||b.in-a.in);
  }
 }
 return {nodes:[...nodes.values()],edges,columns,total:columns[0].reduce((s,n)=>s+n.out,0)};
}
function highlight(id){
 selectedNode=id;
 svg.querySelectorAll('.link').forEach(e=>{const related=!id||e.dataset.source===id||e.dataset.target===id;e.style.opacity=related?(id?'.67':'.34'):'.045';});
 svg.querySelectorAll('.node').forEach(e=>{const n=view.nodes.find(n=>n.id===e.dataset.id);const related=!id||n.id===id||n.linksIn.some(x=>x.source.id===id)||n.linksOut.some(x=>x.target.id===id);e.style.opacity=related?'1':'.3';});
 const n=id?view.nodes.find(n=>n.id===id):null;
 selection.textContent=n?`${n.label} · ${fmt(Math.max(n.in,n.out))} 基因当量 · 仅突出相邻连接`:'悬停查看数值；点击节点突出相邻连接，再次点击或点击空白处取消。';
}
function render(){
 const all=DATA.levels[currentLevel],chosen=filter.value;
 const paths=all.paths.filter(p=>!chosen||p.taxon_key===chosen);
 view=prepare(paths);selectedNode=null;svg.replaceChildren();
 const maxCount=Math.max(1,...view.columns.map(c=>c.length));
 const width=1450,height=Math.max(650,100+maxCount*35),top=90,bottom=40,gap=12,nodeWidth=15;
 const xs=[35,315,825,1130];
 const scale=Math.min(...view.columns.map(c=>(height-top-bottom-gap*Math.max(0,c.length-1))/(view.total||1)));
 svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
 const summary=el('desc');summary.textContent=`${LEVEL_NAMES[currentLevel]}层级，${view.columns[1].length} 个候选分类群；总流量 ${fmt(view.total)} 个分摊基因当量。箭头方向不表示因果。`;
 const headings=['污染物 指标',LEVEL_NAMES[currentLevel]+'分类群','CAZy level2','EC 编号'];
 const subs=['显著候选关联','基因的分类归属','预测家族注释','预测功能编号'];
 for(let s=0;s<4;s++){text(xs[s],32,headings[s],{class:'column-title'});text(xs[s],52,subs[s],{class:'column-sub'});}
 if(!view.edges.length){metrics.textContent='0 个匹配的候选CAZy基因';selection.textContent='没有可显示的连接';text(width/2,height/2,'当前筛选没有可显示的连接',{class:'empty','text-anchor':'middle'});return;}
 for(let s=0;s<4;s++){
  const col=view.columns[s],used=view.total*scale+Math.max(0,col.length-1)*gap;
  let y=top+(height-top-bottom-used)/2;
  for(const n of col){n.x=xs[s];n.y=y;n.h=Math.max(n.in,n.out)*scale;y+=n.h+gap;}
 }
 // Adjacent edges use a single common scale, guaranteeing equal in/out widths.
 for(const n of view.nodes){
  n.linksOut.sort((a,b)=>a.target.y-b.target.y);n.linksIn.sort((a,b)=>a.source.y-b.source.y);
  let out=n.y,inc=n.y;
  for(const e of n.linksOut){e.sy=out;out+=e.value*scale;}
  for(const e of n.linksIn){e.ty=inc;inc+=e.value*scale;}
 }
 const linkLayer=el('g',{'aria-label':'连接'});
 for(const e of view.edges){
  const x0=e.source.x+nodeWidth,x1=e.target.x,c0=x0+(x1-x0)*.47,c1=x1-(x1-x0)*.47,h=e.value*scale;
  const d=`M${x0},${e.sy} C${c0},${e.sy} ${c1},${e.ty} ${x1},${e.ty} L${x1},${e.ty+h} C${c1},${e.ty+h} ${c0},${e.sy+h} ${x0},${e.sy+h} Z`;
  const fill=e.stage===0?(Number(e.rho)>=0?'#c97750':'#527ba6'):'#8da8b4';
  const ribbon=el('path',{d,fill,class:'link','data-source':e.source.id,'data-target':e.target.id,opacity:'.34'},linkLayer);
  let tip=`${e.source.label} → ${e.target.label}\n分摊基因当量：${fmt(e.value)}`;
  if(e.stage===0)tip+=`\n${Number(e.rho)>=0?'正相关':'负相关'}；ρ=${Number(e.rho).toFixed(3)}；原始 p=${safeP(e.p)}`;
  else tip+='\n源文件的预测注释对应，不表示表达或酶活';
  title(ribbon,tip);
  ribbon.addEventListener('mouseenter',()=>{if(!selectedNode)ribbon.style.opacity='.8';});
  ribbon.addEventListener('mouseleave',()=>{if(!selectedNode)ribbon.style.opacity='.34';});
 }
 const nodeLayer=el('g',{'aria-label':'节点'});
 for(const n of view.nodes){
  const g=el('g',{class:'node','data-id':n.id,tabindex:'0',role:'button','aria-label':`${n.label}，${fmt(Math.max(n.in,n.out))} 基因当量`},nodeLayer);
  el('rect',{x:n.x,y:n.y,width:nodeWidth,height:Math.max(.3,n.h),rx:2,fill:nodeColors[n.stage]},g);
  title(g,`${n.label}\n分摊基因当量：${fmt(Math.max(n.in,n.out))}\n点击突出相邻连接`);
  text(n.x+nodeWidth+8,n.y+n.h/2+4,n.label,{class:'node-label'},g);
  const action=()=>highlight(selectedNode===n.id?null:n.id);
  g.addEventListener('click',event=>{event.stopPropagation();action();});
  g.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();action();}});
  g.addEventListener('mouseenter',()=>{if(!selectedNode){highlight(n.id);selectedNode=null;}});
  g.addEventListener('mouseleave',()=>{if(!selectedNode)highlight(null);});
 }
 metrics.innerHTML=`<strong>${view.columns[1].length}</strong> 个${LEVEL_NAMES[currentLevel]} · <strong>${fmt(view.total)}</strong> 基因当量 · ${view.columns[2].length} 个家族显示类别 · ${view.columns[3].length} 个 EC 显示类别`;
 document.getElementById('group-note').textContent='主图合并低权重家族与 EC 为“其他”，但保留其全部权重；“未注释 EC”保留缺失编号的基因份额。筛选单一分类群时仍使用同一套显示类别。';
 selection.textContent='悬停查看数值；点击节点突出相邻连接，再次点击或点击空白处取消。';
}
function populate(){
 filter.replaceChildren();const allOption=document.createElement('option');allOption.value='';allOption.textContent='全部候选分类群';filter.appendChild(allOption);
 for(const taxon of DATA.levels[currentLevel].taxa){const option=document.createElement('option');option.value=taxon.key;option.textContent=taxon.label;filter.appendChild(option);}
 render();
}
for(const tab of document.querySelectorAll('.tab'))tab.addEventListener('click',()=>{currentLevel=tab.dataset.level;for(const t of document.querySelectorAll('.tab'))t.setAttribute('aria-selected',String(t===tab));populate();});
filter.addEventListener('change',render);svg.addEventListener('click',()=>highlight(null));
document.addEventListener('keydown',e=>{if(e.key==='Escape')highlight(null);});
populate();
</script>
</body>
</html>'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    data = json.loads(args.input.read_text(encoding="utf-8"))
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    payload = payload.replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    args.output.write_text(HTML.replace("__DATA__", payload), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
