# NTSR1 多肽设计与筛选

面向小鼠 NTSR1（UniProt O88319）的多肽计算设计流程，结合 PepMLM 条件生成、Boltz 结构预测和 AfDesign 局部优化，完成序列预筛、结构评价、亲本与变体比较、多样性优选及口腔酶切风险注释。

本仓库公开完整流程代码、默认配置、一键启动脚本及可视化工作台。输入数据、模型权重和计算结果由使用者单独准备或生成，准备方式见 [输入说明](docs/INPUTS.md)。本项目调用公开预训练模型，不训练或微调模型权重。

## 功能

- **11 阶段流程**：从靶标准备到候选优选及酶切风险注释。
- **一键运行**：读取已有阶段结果，自动执行缺失步骤，逐阶段输出进度与摘要。
- **可视化工作台**：配置 36 项关键参数及阈值，预览流程、启动任务、查看日志与候选结果。
- **运行记录**：保存配置、随机种子、阶段日志；统一入口额外记录运行环境与模型文件哈希。

## 目录

```text
.
├── main.py                 # 全流程运行与标准结果导出
├── run_pipeline.sh         # 复用结果并运行缺失阶段
├── run_*.sh                # 分阶段运行入口
├── workbench.py            # 工作台服务
├── pipeline/               # 算法及调度代码
├── config/                 # 默认参数与随机种子
├── environment/            # 环境及依赖配置
├── templates/              # 工作台页面
├── workbench_static/       # 页面样式、脚本与 Logo
├── docs/                   # 输入、模型及使用说明
├── tests/                  # 不依赖研究数据的工作台测试
└── LICENSE                 # MIT 许可证
```

所有命令均在仓库根目录执行。

## 运行环境

工作台及执行计划预览使用 Python 3.10 或以上。模型计算使用 Linux x86_64、兼容 CUDA 的 NVIDIA GPU 和三个 Conda 环境：

|环境|Python|主要组件|
|---|---|---|
|pepmlm|3.10|PepMLM-650M、PyTorch、Transformers|
|afdesign|3.10|ColabDesign v1.1.3、JAX CUDA 12|
|boltz|3.11|Boltz 2.2.1、PyTorch 2.6.0/cu124|

```bash
conda env create -f environment/pepmlm.yml
conda env create -f environment/afdesign.yml
conda env create -f environment/boltz.yml
```

原计算环境参考：Ubuntu 22.04、NVIDIA A800 80 GB、64 GB 主存及至少 100 GB 可用磁盘。CUDA 12.4 对应驱动要求见 [NVIDIA 发布说明](https://docs.nvidia.com/cuda/archive/12.4.0/cuda-toolkit-release-notes/index.html)。环境文件列出依赖版本或允许范围；完整模型计算耗时取决于设备、生成规模和预测数量。

## 快速开始

### 预览流程

无需模型权重或输入数据即可生成执行计划：

```bash
python main.py --mode full --plan-only
```

### 运行筛选

按 [输入准备说明](docs/INPUTS.md) 准备参考包、受体 MSA、酶切统计和模型权重，然后执行：

```bash
bash run_pipeline.sh
```

脚本检查 `work/` 中各阶段输出，复用已有结果并计算缺失步骤。仅展示已有结果或按新参数重新计算：

```bash
bash run_pipeline.sh --summary-only
bash run_pipeline.sh --work-dir work/new_run
```

统一入口可从头运行，并自动导出结果：

```bash
python main.py --mode full --output-dir work/reproduce
```

该方式将配置、日志、环境记录和筛选结果写入 `work/reproduce/`，正式候选表位于 `work/reproduce/results/results.csv`。各阶段也可通过 `pipeline/` 中的脚本及对应 `run_*.sh` 单独执行。

### 启动工作台

```bash
python -m pip install -r environment/requirements-workbench.txt
python workbench.py --port 8081
```

使用 Cursor Remote SSH 时，转发端口 `8081` 后在浏览器打开。未产生筛选结果时，页面显示空结果；可以先配置参数并预览流程。计算任务需要先准备模型环境和输入文件。

工作台依次提供运行日志、流程总览和候选结果，支持检索候选及下载 CSV、PDB 和日志。完整操作说明见 [工作台教程](docs/WORKBENCH.md)。

## 筛选流程

|阶段|内容|
|---|---|
|01|准备靶标序列、参考结构与口袋掩码|
|02|PepMLM 条件生成及完整序列去重|
|03|理化过滤与 PPL 计算|
|04|肽长分层、PPL 排序与聚类去冗余|
|05–06|初始复合物结构任务、Boltz 预测与结构筛选|
|07|AfDesign 局部序列优化|
|08–09|变体独立结构预测及亲本配对复核|
|10|结合改进、退化保护与序列多样性优选候选|
|11|口腔酶切风险注释，保持原候选排序|

默认设置生成 100,000 条、肽长 12–33 aa、最终候选数 100；这些数量可在配置或工作台中修改。当前流程针对小鼠 NTSR1，替换靶标时需同步调整参考输入和相关配置。

## 输出与说明

阶段结果写入 `work/` 下的对应目录。`main.py --mode export` 将已生成的阶段 10/11 结果合并为 `results/results.csv`，同时保留候选结构和预测证据。结果表包含候选编号、序列、排名、结构指标、酶切风险、模型与运行版本等字段。

结构置信度和序列评分用于计算筛选；酶切风险是候选集合内的相对注释，不表示实测切割概率或半衰期。MD 与湿实验不属于本代码流程。详见 [模型与输入说明](docs/DATA_MODEL_CARD.md) 及 [酶切风险说明](docs/11_口腔酶切风险初筛说明.md)。

## 许可

自主开发的筛选流程与可视化工作台代码采用 [MIT 许可证](LICENSE)。第三方模型、权重、数据及图片等资源遵循各自许可，来源见 [模型与输入说明](docs/DATA_MODEL_CARD.md)。
