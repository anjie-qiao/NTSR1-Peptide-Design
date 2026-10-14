# 输入准备

代码版提供脚本与配置。运行模型前，需单独准备以下输入；候选表、结构、运行日志和中间结果由流程生成。

## 1. 文件位置

|输入|默认位置|用途|
|---|---|---|
|NTSR1 参考包|`models/mouse_NTSR1_6Z8N_6ZIN_pack.zip`|阶段 01 的序列、结构映射和口袋掩码|
|受体 MSA|`models/ntsr1_msa/ntsr1.csv`|两轮 Boltz 预测的受体输入|
|MEROPS 统计快照|`data/oral_proteolysis/merops_specificity_20260929.json`|阶段 11 的位点偏好统计|
|PepMLM 配置、tokenizer 和权重|`models/PepMLM-650M/`|候选生成与 PPL 计算|
|AlphaFold Multimer v3 参数|`models/afdesign/params/params_model_1_multimer_v3.npz`|AfDesign 优化|
|Boltz confidence/affinity 权重和 CCD|`models/boltz_cache/`|结构预测|

前三项属于项目运行输入，当前代码发布不附带这些文件。参考包是经过序列映射、结构加工和口袋整理的组合输入，不能以单个 RCSB 原始结构文件直接替代。准备自己的输入时，应按以下结构组织，并在 `config/` 中更新路径。

## 2. NTSR1 参考包

序列来源为 [UniProt O88319](https://www.uniprot.org/uniprotkb/O88319/entry)，参考结构来源为 [RCSB 6Z8N](https://www.rcsb.org/structure/6Z8N) 和 [6ZIN](https://www.rcsb.org/structure/6ZIN)。ZIP 内部文件名应符合阶段 01 的 `MODEL_SPECS`：

```text
SHA256SUMS
source_info/
  O88319_mouse_NTSR1.fasta
  mouse_NTSR1_6Z8N_missing_residues.csv
  mouse_NTSR1_6ZIN_missing_residues.csv
models/
  mouse_NTSR1_6Z8N.pdb
  mouse_NTSR1_6Z8N.cif
  mouse_NTSR1_6ZIN.pdb
  mouse_NTSR1_6ZIN.cif
ligands/
  6Z8N_SR9_original_coords.pdb
  6ZIN_Q6Q_original_coords.pdb
pockets/
  mouse_NTSR1_6Z8N_pocket5A.csv
  mouse_NTSR1_6ZIN_pocket5A.csv
```

- FASTA 标题包含 `O88319`，序列为标准氨基酸；默认长度为 424，见阶段 01 配置。
- 结构使用受体链 A，残基编号和序列应映射至小鼠 NTSR1 的标准编号。
- 口袋 CSV 至少包含 `one_letter`、`mouse_residue_number`。
- 缺失残基 CSV 包含 `missing_mouse_residue_start`、`missing_mouse_residue_end`、`count`。
- `SHA256SUMS` 每行记录文件的 SHA256 和 ZIP 内相对路径。阶段 01 会验证这些文件。

## 3. 受体 MSA

准备符合 Boltz 输入格式的 CSV，包含 `key`、`sequence` 两列，并与当前受体序列一致。其位置由阶段 05 和 08 的 `receptor_msa` 指定。肽链默认使用 `msa: empty`；本流程不调用在线 MSA 服务。

## 4. 酶切统计

从 [MEROPS](https://www.ebi.ac.uk/merops/) 的底物位点统计整理快照，保存数据来源、获取日期及适用范围。默认脚本使用五个 profile：

|profile|MEROPS 编号|
|---|---|
|neutrophil_elastase|S01.131|
|cathepsin_B|C01.060|
|cathepsin_L|C01.032|
|MMP8|M10.002|
|MMP9|M10.004|

JSON 必须包含 `positions`、`scope`、`entries`。`positions` 顺序为 `P4, P3, P2, P1, P1prime, P2prime, P3prime, P4prime`；`entries` 按 MEROPS 编号索引，每个 profile 包含 `source` 和 `counts`。

`counts` 为 20 种标准氨基酸到 8 位计数数组的映射，值为非负整数，每个位点至少存在一个非零观测。可额外记录 `schema_version`、`retrieved_date`、profile 名称及观测总数。当前发布不包含原统计快照或替代统计值。

## 5. 模型权重

在具备网络访问的环境中下载后，也可将模型文件复制到计算服务器：

```bash
conda run -n pepmlm hf download ChatterjeeLab/PepMLM-650M \
  --local-dir models/PepMLM-650M

mkdir -p models/afdesign/params models/boltz_cache
curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar \
  | tar -x -C models/afdesign/params
wget -c https://model-gateway.boltz.bio/boltz2_conf.ckpt -P models/boltz_cache
wget -c https://model-gateway.boltz.bio/boltz2_aff.ckpt -P models/boltz_cache
wget -c https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar -P models/boltz_cache
tar -xf models/boltz_cache/mols.tar -C models/boltz_cache
```

Boltz 缓存需包含两个 `.ckpt`、`mols.tar` 和解压后的 `mols/`。affinity 权重用于模型初始化，本流程不输出亲和力预测。来源与许可见 [模型说明](DATA_MODEL_CARD.md)。
