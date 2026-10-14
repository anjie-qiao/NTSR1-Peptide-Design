'use strict';
(() => {
  const root = document.getElementById('ntsr1-demo');
  const $ = id => document.getElementById(id);
  const csrf = document.querySelector('meta[name="workbench-csrf"]').content;
  let config = null, state = null, source = 'historical', view = 'flow', selected = 4, workspace = 'results';
  let values = {}, rows = [], page = 0, chosen = null, polling = false, submitting = false;
  const labels = {queued:'等待启动', running:'正在运行', planned:'流程预览已生成', complete:'筛选完成', failed:'运行失败', interrupted:'进程中断'};
  const fmt = value => value == null ? '—' : Number(value).toLocaleString('zh-CN');
  const number = (value, digits=4) => value === '' || value == null ? '—' : Number(value).toFixed(digits);
  const download = (kind, candidate) => '/api/download/' + kind + '?' + new URLSearchParams({source, ...(candidate ? {candidate} : {})});

  async function api(path, data) {
    const response = await fetch(path, data === undefined ? {} : {method:'POST',headers:{'Content-Type':'application/json','X-Workbench-CSRF':csrf},body:JSON.stringify(data)});
    const result = await response.json();
    if (!response.ok) { if (response.status === 401) $('login').showModal(); throw new Error(result.error || '操作未完成，请查看服务器日志。'); }
    return result;
  }
  function notice(message, error=false) { $('notice').textContent=message; $('notice').hidden=!message; $('notice').classList.toggle('error',error); }
  function dirty() { return config && Object.keys(values).some(key => values[key] !== config.values[key]); }
  function updateDirty() {
    const busy=submitting||state?.tasks.some(t=>['queued','running'].includes(t.status));
    const changed=dirty();
    $('config-state').textContent=busy?'任务运行中':changed?'有未保存的修改':'参数已载入';
    $('save').disabled=busy||!config||!changed;
  }
  function bindInput(input, key) {
    input.addEventListener('input',()=>{
      values[key]=input.value === '' ? null : Number(input.value);
      updateDirty();
    });
  }
  function populateConfig(data) {
    config=data; values={...data.values};
    $('advanced-fields').replaceChildren();
    $('parameter-index').replaceChildren();
    const groups = new Map();
    for (const f of data.fields) {
      if (!groups.has(f.group)) {
        const section=document.createElement('section'), heading=document.createElement('div'), number=document.createElement('span'), title=document.createElement('h3'), grid=document.createElement('div'), link=document.createElement('a');
        const parts=f.group.split(' · ');
        section.className='parameter-group';section.id='parameter-group-'+groups.size;section.dataset.stage=f.group==='基本参数'?'basic':f.stage;
        heading.className='parameter-heading';number.className='stage-index';number.textContent=parts.length===2?parts[0]:'01';title.textContent=parts.at(-1);heading.append(number,title);
        grid.className='parameter-grid';section.append(heading,grid);$('advanced-fields').append(section);groups.set(f.group,grid);
        link.href='#'+section.id;link.textContent=parts.at(-1);$('parameter-index').append(link);
      }
      const label=document.createElement('label'), title=document.createElement('span'), input=document.createElement('input');
      label.className='field';title.textContent=f.label;
      if(f.nullable){const optional=document.createElement('span');optional.className='optional-field';optional.textContent='可选';title.append(optional);}
      input.id='parameter-'+f.key;input.type='number';input.min=f.min;input.max=f.max;input.step=f.step;input.dataset.param=f.key;input.value=values[f.key] ?? '';input.required=!f.nullable;input.placeholder=f.nullable?'不启用':'';bindInput(input,f.key);label.append(title,input);groups.get(f.group).append(label);
    }
    updateDirty();
  }
  async function loadConfig() { populateConfig(await api('/api/config')); }
  function setWorkspace(next) {
    workspace=next;
    $('workspace-results').hidden=next!=='results';$('workspace-parameters').hidden=next!=='parameters';
    $('runs-open').setAttribute('aria-pressed',String(next==='results'));$('advanced-open').setAttribute('aria-pressed',String(next==='parameters'));
  }
  function openParameters(stage) {
    setWorkspace('parameters');
    const target=stage?$('advanced-fields').querySelector('[data-stage="'+stage+'"]'):window.innerWidth<=680?$('workspace-parameters'):null;
    target?.scrollIntoView({block:'start'});
  }
  $('advanced-open').addEventListener('click',()=>openParameters());
  $('runs-open').addEventListener('click',()=>setWorkspace('results'));
  $('save').addEventListener('click',async()=>{
    $('save').disabled=true;
    try {
      const changed=Object.fromEntries(Object.entries(values).filter(([key,value])=>value!==config.values[key]));
      if(!Object.keys(changed).length){notice('参数未发生变化。');return;}
      const result=await api('/api/config',{revision:config.revision,values:changed});
      populateConfig({...config,...result});notice('参数已保存。');await refresh();
    } catch(e){notice(e.message,true);} finally{updateDirty();}
  });
  async function launchTask(mode) {
    if(submitting)return;
    if(dirty()){openParameters();notice('请先保存参数。',true);return;}
    submitting=true;updateButtons();
    try {
      const job=await api('/api/tasks',{name:$('task-name').value,mode,revision:config.revision});
      source=job.id;selected=0;rows=[];chosen=null;page=0;notice('任务已启动。');setWorkspace('results');setView('logs');await refresh();
    }catch(e){notice(e.message,true);}finally{submitting=false;if(state)updateButtons();}
  }
  $('preview-run').addEventListener('click',()=>launchTask('plan'));
  $('run').addEventListener('click',()=>launchTask('full'));
  function setView(next) {
    view=next;
    for(const name of ['flow','results','logs']){$(name).hidden=name!==next;$(name+'-tab').setAttribute('aria-selected',String(name===next));}
    if(next==='results')loadResults().catch(e=>notice(e.message,true));
    if(next==='logs')loadLog().catch(e=>notice(e.message,true));
  }
  for(const name of ['flow','results','logs'])$(name+'-tab').addEventListener('click',()=>setView(name));
  $('source').addEventListener('change',async()=>{source=$('source').value;rows=[];chosen=null;page=0;notice('');await refresh();if(view==='results')await loadResults();if(view==='logs')await loadLog();});
  function updateButtons() {
    const busy=submitting||state.tasks.some(t=>['queued','running'].includes(t.status));
    $('preview-run').disabled=busy||!config;
    $('run').disabled=busy||!config||!state.environment.full_available;
    root.querySelectorAll('[data-param]').forEach(input=>input.disabled=busy);
    updateDirty();
  }
  function detailRows(items) {
    const list=document.createElement('ul');list.className='rules';
    for(const [label,value] of items){const li=document.createElement('li'), span=document.createElement('span'), strong=document.createElement('span');span.className='small';span.textContent=label;strong.className='rule-value';strong.textContent=value;li.append(span,strong);list.append(li);}return list;
  }
  function detailsFor(i) {
    const p=state.parameters, s=state.stages[i];
    return [
      [['靶标','小鼠 NTSR1 · O88319'],['参考结构','6Z8N / 6ZIN'],['输入','项目预置参考包']],
      [['生成规模',fmt(p.total_candidates)+' 条'],['肽长',p.length_min+'–'+p.length_max+' aa'],['采样','四种配置 · 固定随机种子']],
      [['理化通过',fmt(s.count)+' 条'],['进入 PPL 计算',fmt(s.counts.ppl_computed)+' 条'],['净电荷范围',p.qc_min_net_charge+' 至 '+p.qc_max_net_charge],['每肽长 PPL 数量上限',fmt(p.ppl_limit)]],
      [['PPL 保留比例',(p.ppl_keep*100).toFixed(0)+'%'],['聚类一致性阈值',(p.cluster_identity*100).toFixed(0)+'%'],['核心库',fmt(s.count)+' 条']],
      [['平均 pair ipTM','≥ '+p.structure_min_pair_iptm_mean],['平均界面 PAE','≤ '+p.structure_max_interface_pae_mean+' Å'],['肽接触比例','≥ '+p.structure_min_peptide_contact_fraction],['严重原子冲突','≤ '+p.structure_max_severe_clash_atom_pairs+' 对']],
      [['亲本家族',fmt(s.counts.input_parents)+' 个'],['累计优化','一位 / 二位 / 三位替换'],['生成变体',fmt(s.counts.generated_variants)+' 个'],['序列通过变体',fmt(s.count)+' 个'],['比较方式','逐亲本配对比较']],
      [['目标候选数',fmt(p.final_count)+' 条'],['已产生候选',fmt(s.count)+' 条'],['家族优选','合格变体与亲本配对比较'],['最终一致性阈值',String(p.final_identity)]],
      [['候选注释数量',fmt(s.count)+' 条'],['风险分层','相对低 / 中 / 高'],['结果','候选表 / 酶切注释表 / PDB']]
    ][i];
  }
  function showStage(i) {
    selected=i;root.querySelectorAll('[data-stage-index]').forEach(b=>{const active=Number(b.dataset.stageIndex)===i;b.classList.toggle('selected',active);b.setAttribute('aria-pressed',String(active));});
    const s=state.stages[i], detail=$('stage-detail');detail.replaceChildren();
    const eyebrow=document.createElement('div'), title=document.createElement('h3'), note=document.createElement('div');eyebrow.className='eyebrow';eyebrow.textContent='筛选依据';title.textContent=s.title;note.className='small';note.textContent=source==='historical'?'当前配置':'任务参数';detail.append(eyebrow,title,note,detailRows(detailsFor(i)));
    const stage=[1,2,3,4,6,7,10,11][i];
    if([3,4,6,10].includes(stage)){const edit=document.createElement('button');edit.type='button';edit.textContent='配置此阶段参数';edit.addEventListener('click',()=>openParameters(stage));detail.append(edit);}
  }
  function renderStages() {
    $('stages').replaceChildren();
    const job=state.task;
    for(const [i,s] of state.stages.entries()) {
      const button=document.createElement('button');button.type='button';button.className='stage';button.dataset.stageIndex=i;
      const index=document.createElement('span'), text=document.createElement('span'), title=document.createElement('span'), sub=document.createElement('span'), amount=document.createElement('span'), status=document.createElement('span');
      index.className='stage-index';index.textContent=String(i+1).padStart(2,'0');title.className='stage-name';title.textContent=s.title;sub.className='stage-sub';sub.textContent=s.subtitle;text.append(title,sub);
      amount.className='stage-value';amount.textContent=fmt(s.count)+(i===0&&s.count!=null?' aa':'');status.className='stage-count-state';
      const startStage=[1,2,3,4,5,7,10,11][i];
      const active=job?.status==='running'&&job.active_stage>=startStage&&job.active_stage<=s.stage;
      status.textContent=active?'运行中':s.count!=null?'已完成':'待运行';status.hidden=source==='historical'&&s.count!=null;button.classList.toggle('active',active);amount.append(status);button.append(index,text,amount);button.addEventListener('click',()=>showStage(i));$('stages').append(button);
    }
    showStage(selected);
  }
  async function refresh() {
    state=await api('/api/state?'+new URLSearchParams({source}));
    const hasResults=state.result_count>0||state.stages.some(s=>s.count!=null), allResults=state.stages.every(s=>s.count!=null);
    const selector=$('source');selector.replaceChildren(new Option(hasResults?'NTSR1 · 已有结果':'NTSR1 · 待运行','historical'));
    for(const t of state.tasks)selector.append(new Option(t.name+' · '+labels[t.status],t.id));selector.value=source;
    const job=state.task; $('status-badge').textContent=job?labels[job.status]:allResults?'筛选完成':hasResults?'已有阶段结果':'待运行';
    $('status').textContent=job?(job.mode==='plan'?(job.status==='planned'?'流程预览已生成':labels[job.status]+' · 流程预览'):(labels[job.status]+(job.active_stage?' · 阶段 '+String(job.active_stage).padStart(2,'0')+(job.phase==='boltz'?' · Boltz 预测':''):''))):allResults?'11 / 11 阶段完成':hasResults?'阶段结果已载入':'尚无运行结果';
    if(job?.error)notice(job.error,true);
    const completed=source==='historical'?(allResults?11:0):job?.status==='complete'?11:job?.mode==='plan'?0:Math.max(0,(job?.active_stage||1)-1);
    $('progress').style.width=completed/11*100+'%';root.querySelector('.track').setAttribute('aria-valuenow',completed);
    $('platform').textContent=state.environment.platform;
    $('environment-detail').textContent=[...state.environment.reasons,...state.environment.missing,state.environment.note].join('；');
    updateButtons();renderStages();
    $('download-plan').hidden=job?.status!=='planned';$('download-plan').href=download('plan');$('download-log').hidden=!job;$('download-log').href=download('log');
    if(workspace==='results'&&view==='logs')await loadLog();
    if(workspace==='results'&&view==='results')await loadResults();
  }
  async function loadLog() {
    if(source==='historical'){$('log-text').textContent='选择一个运行任务查看日志。';$('log-note').textContent='';return;}
    const result=await api('/api/log?'+new URLSearchParams({source}));$('log-text').textContent=result.text;$('log-note').textContent=result.truncated?'显示最近 64 KB；可下载完整日志。':'';
  }
  async function loadResults() {
    const result=await api('/api/results?'+new URLSearchParams({source}));rows=result.rows;
    $('results-title').textContent='候选结果 · '+fmt(rows.length)+' 条';$('download-results').hidden=!rows.length;$('download-results').href=download('results');renderResults();
  }
  function filteredRows(){const query=$('search').value.trim().toUpperCase();return rows.filter(r=>r.candidate_id.toUpperCase().includes(query)||r.sequence.toUpperCase().includes(query));}
  function renderResults() {
    const filtered=filteredRows();page=Math.min(page,Math.max(0,Math.ceil(filtered.length/10)-1));$('candidate-rows').replaceChildren();
    for(const r of filtered.slice(page*10,page*10+10)) {
      const tr=document.createElement('tr'), rank=document.createElement('td'), identifier=document.createElement('td'), button=document.createElement('button'), type=document.createElement('span');rank.textContent=r.final_rank;button.className='candidate-btn';button.type='button';button.textContent=r.candidate_id;button.addEventListener('click',()=>showCandidate(r));type.className='candidate-type';type.textContent=r.source_type==='afdesign'?'优化变体':'原始亲本';identifier.append(button,type);tr.append(rank,identifier);
      for(const [key,digits] of [['pair_iptm_mean',4],['interface_pae_mean',2]]){const td=document.createElement('td');td.className='num';td.textContent=number(r[key],digits);tr.append(td);}
      const risk=document.createElement('td'),level=document.createElement('span'),riskSource=document.createElement('td');
      const riskValue=r.oral_cleavage_risk_level||'';
      level.className='risk-level'+(riskValue?(/低|low/i.test(riskValue)?' low':/高|high/i.test(riskValue)?' high':' medium'):'');level.textContent=riskValue||'—';risk.append(level);
      if(r.oral_cleavage_risk_percentile!=null&&r.oral_cleavage_risk_percentile!==''){const percentile=document.createElement('span');percentile.className='risk-percentile';percentile.textContent='百分位 '+number(r.oral_cleavage_risk_percentile,2);risk.append(percentile);}
      riskSource.className='risk-source';riskSource.textContent=r.main_risk_source||'—';tr.append(risk,riskSource);$('candidate-rows').append(tr);
    }
    if(!filtered.length){const tr=document.createElement('tr'),td=document.createElement('td');td.colSpan=6;td.textContent=rows.length?'没有匹配候选。':'此任务尚未生成候选结果。';tr.append(td);$('candidate-rows').append(tr);}
    $('page-label').textContent=filtered.length?`${page+1} / ${Math.ceil(filtered.length/10)} 页 · ${filtered.length} 条`:'0 条';$('previous').disabled=page===0;$('next').disabled=(page+1)*10>=filtered.length;
    if(chosen){const candidate=rows.find(r=>r.candidate_id===chosen);if(candidate)showCandidate(candidate);else $('candidate-detail').hidden=true;}else if(rows.length)showCandidate(rows[0]);else $('candidate-detail').hidden=true;
  }
  function showCandidate(r) {
    chosen=r.candidate_id;const panel=$('candidate-detail');panel.hidden=false;panel.replaceChildren();
    const title=document.createElement('h3'),seq=document.createElement('div'),note=document.createElement('div');title.textContent=r.candidate_id;title.style.overflowWrap='anywhere';seq.className='sequence';seq.textContent=r.sequence;note.className='small';note.textContent=`排名 ${r.final_rank} · ${r.length} aa · ${r.source_type==='afdesign'?'优化变体':'原始亲本'}`;
    panel.append(title,seq,note,detailRows([['平均肽 pLDDT',number(r.peptide_mean_plddt)],['肽接触比例',r.peptide_contact_fraction?(Number(r.peptide_contact_fraction)*100).toFixed(1)+'%':'—'],['相对酶切风险',r.oral_cleavage_risk_level||'—'],['酶切风险百分位',number(r.oral_cleavage_risk_percentile,2)],['主要风险来源',r.main_risk_source||'—']]));
    if(r.structure_pdb){const link=document.createElement('a');link.className='download-action';link.textContent='下载对应 PDB 结构';link.href=download('structure',r.candidate_id);panel.append(link);}
  }
  $('search').addEventListener('input',()=>{page=0;renderResults();});$('previous').addEventListener('click',()=>{page--;renderResults();});$('next').addEventListener('click',()=>{page++;renderResults();});
  async function initialise(){await loadConfig();await refresh();}
  $('login-button').addEventListener('click',async()=>{try{await api('/api/login',{token:$('access-token').value});$('access-token').value='';$('login').close();await initialise();}catch(e){$('login-error').textContent=e.message;}});
  if(root.dataset.needsLogin==='true')$('login').showModal();else initialise().catch(e=>notice(e.message,true));
  setInterval(async()=>{if(polling||!state||$('login').open)return;polling=true;try{await refresh();}catch(e){notice(e.message,true);}finally{polling=false;}},2500);
})();
