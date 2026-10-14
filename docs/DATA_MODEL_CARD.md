# 模型与输入说明

本项目使用公开预训练模型，不训练或微调模型权重；AfDesign 优化候选序列。当前发布范围为源码、配置和使用说明，输入数据及模型权重按 [输入准备说明](INPUTS.md) 单独准备。

## 模型

|模型|来源与版本|主要设置|许可|
|---|---|---|---|
|PepMLM-650M|[ChatterjeeLab/PepMLM-650M](https://huggingface.co/ChatterjeeLab/PepMLM-650M)，约 6.5 亿参数|受体序列条件生成、逐位掩码 PPL|模型卡标注 MIT|
|Boltz-2|[Boltz 2.2.1](https://github.com/jwohlwend/boltz/tree/v2.2.1)，PyTorch 2.6.0/cu124|3 recycles、200 sampling steps、1 diffusion sample、max_msa_seqs=512、use_potentials、full PAE、no_kernels|[代码及权重 MIT](https://github.com/jwohlwend/boltz#license)|
|AfDesign|[ColabDesign v1.1.3](https://github.com/sokrypton/ColabDesign/tree/v1.1.3)，AlphaFold 参数 2022-12-06，model_1_multimer_v3|1 recycle、1 model、dropout=false；三轮搜索，每轮 3 proposals|ColabDesign [Beer-Ware](https://github.com/sokrypton/ColabDesign/blob/v1.1.3/LICENSE.txt)；AlphaFold 代码 Apache 2.0、参数 [CC BY 4.0](https://github.com/google-deepmind/alphafold#license-and-disclaimer)|

参考项目：[PepMLM](https://github.com/programmablebio/pepmlm)。代码调用模型接口，默认 `validate_upstream_count=0`，不附上游可选生成脚本。

## 输入来源

|输入|来源|用途与许可|
|---|---|---|
|小鼠 NTSR1 序列|[UniProt O88319](https://www.uniprot.org/uniprotkb/O88319/entry)|统一受体序列输入；[CC BY 4.0](https://www.uniprot.org/help/license)|
|结构参考|[RCSB 6Z8N](https://www.rcsb.org/structure/6Z8N)、[6ZIN](https://www.rcsb.org/structure/6ZIN)|需映射并整理为阶段 01 参考包；遵循 [RCSB 政策](https://www.rcsb.org/pages/policies)|
|受体 MSA|使用者准备与目标序列一致的 MSA|Boltz 受体输入，按所用数据来源的条款使用|
|酶切位点统计|[MEROPS](https://www.ebi.ac.uk/merops/)|阶段 11 的五个 profile；[EMBL-EBI 使用条款](https://www.ebi.ac.uk/about/terms-of-use/)|
|酶切基序规则|参考文献链接记录在 `config/11_oral_proteolysis_filter.json`|RgpA/RgpB、Kgp、DPP-IV 的序列基序注释|

## 默认参数与去重

生成 seed=20260810、AfDesign seed=20260817；自动运行入口的 Boltz seed=42。生成 profile、理化阈值和结构门槛见 `config/`。

阶段 02 按完整序列去重；阶段 04 按肽长分层，并按 PPL 排序和等长位置一致性聚类；阶段 10 每个亲本家族保留一个代表，再按序列多样性优选。默认阶段 04 一致性阈值为 0.65，阶段 10 为 0.70。

本项目不开展模型训练，因此没有项目内训练、验证及测试集划分。候选去重针对本次设计集合；优化变体由 Boltz 独立复核，并与对应亲本配对比较。

`main.py` 完整运行时记录模型文件 SHA256、已安装环境版本和阶段耗时，写入运行目录的 `runtime_environment.json` 与日志。

## 适用范围

当前输入与结构编号针对小鼠 NTSR1。PepMLM 使用受体序列；PPL、结构置信度和接触指标用于计算筛选，不直接代表亲和力、特异性或受体功能。

结构预测不强制使用 6Z8N/6ZIN 模板；口袋接触比例参与指标计算与排序。膜环境、胞外可达性、糖基化和脱靶作用需要进一步评价。

酶切注释基于位点统计和序列基序，给出候选集合内的相对风险，不改变阶段 10 排名；不包含 pH、温度、端基修饰及结构可及性建模。MD 和湿实验不在本代码流程中。
