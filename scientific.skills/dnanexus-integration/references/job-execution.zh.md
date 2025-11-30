<!-- 此文件由机器翻译自 job-execution.md -->

# DNAnexus 作业执行和工作流程

## 概述

作业是 DNAnexus 上的基本执行单元。当小程序或应用程序运行时，会在具有持续 API 访问的隔离 Linux 环境中的工作节点上创建并执行作业。

## 工作类型

### 起源工作
最初由用户或自动化系统创建。

### 大师职位
直接启动可执行文件（应用程序/小程序）的结果。

### 儿童工作
由用于并行处理或子工作流的父作业生成。

## 运行作业

### 运行小程序

**基本执行**：
```python
import dxpy

# Run an applet
job = dxpy.DXApplet("applet-xxxx").run({
    "input1": {"$dnanexus_link": "file-yyyy"},
    "input2": "parameter_value"
})

print(f"Job ID: {job.get_id()}")
```

**使用命令行**：
<<<代码块_1>>>

### 运行应用程序

<<<代码块_2>>>

### 指定执行参数

<<<代码块_3>>>

## 作业监控

### 检查作业状态

<<<代码块_4>>>

**使用命令行**：
<<<代码块_5>>>

### 等待作业完成

<<<代码块_6>>>

### 获取作业输出

```python
job = dxpy.DXJob("job-xxxx")

# Wait for completion
job.wait_on_done()

# Get outputs
output = job.describe()["output"]
output_file_id = output["result_file"]["$dnanexus_link"]

# Download result
dxpy.download_dxfile(output_file_id, "result.txt")
```

### 作业输出参考

在作业输出完成之前创建对作业输出的引用：

```python
# Launch first job
job1 = dxpy.DXApplet("applet-1").run({"input": "..."})

# Launch second job using output reference
job2 = dxpy.DXApplet("applet-2").run({
    "input": dxpy.dxlink(job1.get_output_ref("output_name"))
})
```

## 作业日志

### 查看日志

**命令行**：
```bash
dx watch job-xxxx --get-streams
```

**以编程方式**：
```python
import sys

# Get job logs
job = dxpy.DXJob("job-xxxx")
log = dxpy.api.job_get_log(job.get_id())

for log_entry in log["loglines"]:
    print(log_entry)
```

## 并行执行

### 创建子作业

```python
@dxpy.entry_point('main')
def main(input_files):
    # Create subjobs for parallel processing
    subjobs = []

    for input_file in input_files:
        subjob = dxpy.new_dxjob(
            fn_input={"file": input_file},
            fn_name="process_file"
        )
        subjobs.append(subjob)

    # Collect results
    results = []
    for subjob in subjobs:
        result = subjob.get_output_ref("processed_file")
        results.append(result)

    return {"all_results": results}

@dxpy.entry_point('process_file')
def process_file(file):
    # Process single file
    # ...
    return {"processed_file": output_file}
```

### 分散-聚集模式

```python
# Scatter: Process items in parallel
scatter_jobs = []
for item in items:
    job = dxpy.new_dxjob(
        fn_input={"item": item},
        fn_name="process_item"
    )
    scatter_jobs.append(job)

# Gather: Combine results
gather_job = dxpy.new_dxjob(
    fn_input={
        "results": [job.get_output_ref("result") for job in scatter_jobs]
    },
    fn_name="combine_results"
)
```

## 工作流程

工作流将多个应用程序/小程序组合成多步骤管道。

### 创建工作流程

```python
# Create workflow
workflow = dxpy.new_dxworkflow(
    name="My Analysis Pipeline",
    project="project-xxxx"
)

# Add stages
stage1 = workflow.add_stage(
    dxpy.DXApplet("applet-1"),
    name="Quality Control",
    folder="/qc"
)

stage2 = workflow.add_stage(
    dxpy.DXApplet("applet-2"),
    name="Alignment",
    folder="/alignment"
)

# Connect stages
stage2.set_input("reads", stage1.get_output_ref("filtered_reads"))

# Close workflow
workflow.close()
```

### 运行工作流程

```python
# Run workflow
analysis = workflow.run({
    "stage-xxxx.input1": {"$dnanexus_link": "file-yyyy"}
})

# Monitor analysis (collection of jobs)
analysis.wait_on_done()

# Get workflow outputs
outputs = analysis.describe()["output"]
```

**使用命令行**：
```bash
dx run workflow-xxxx -i stage-1.input=file-yyyy
```

## 作业权限和上下文

### 工作区上下文

作业在具有克隆输入数据的工作区项目中运行：
- 作业需要 `CONTRIBUTE` 工作区权限
- 作业需要 `VIEW` 访问源项目
- 所有费用均计入原始项目

### 数据要求

作业无法开始，直到：
1. 所有输入数据对象均处于`closed`状态
2. 所需权限已具备
3. 资源分配

在工作区清理之前，输出对象必须达到 `closed` 状态。

## 作业生命周期

```
Created → Waiting on Input → Runnable → Running → Done/Failed
```

**国家**：
- `idle`：作业已创建但尚未排队
- `waiting_on_input`：等待输入数据对象关闭
- `runnable`：准备运行，等待资源
- `running`：当前正在执行
- `done`：成功完成
- `failed`：执行失败
- `terminated`：手动停止

## 错误处理

### 工作失败

```python
job = dxpy.DXJob("job-xxxx")
job.wait_on_done()

desc = job.describe()
if desc["state"] == "failed":
    print(f"Job failed: {desc.get('failureReason', 'Unknown')}")
    print(f"Failure message: {desc.get('failureMessage', '')}")
```

### 重试失败的作业

```python
# Rerun failed job
new_job = dxpy.DXApplet(desc["applet"]).run(
    desc["originalInput"],
    project=desc["project"]
)
```

### 终止作业

```python
# Stop a running job
job = dxpy.DXJob("job-xxxx")
job.terminate()
```

**使用命令行**：
```bash
dx terminate job-xxxx
```

## 资源管理

### 实例类型

指定计算资源：

```python
# Run with specific instance type
job = dxpy.DXApplet("applet-xxxx").run(
    {"input": "..."},
    instance_type="mem3_ssd1_v2_x8"  # 8 cores, high memory, SSD
)
```

常见实例类型：
- `mem1_ssd1_v2_x4` - 4 核，标准内存
- `mem2_ssd1_v2_x8` - 8 核，高内存
- `mem3_ssd1_v2_x16` - 16 核，非常高的内存
- `mem1_ssd1_v2_x36` - 36 个内核用于并行工作负载

### 超时设置

设置最大执行时间：

```python
job = dxpy.DXApplet("applet-xxxx").run(
    {"input": "..."},
    timeout="24h"  # Maximum runtime
)
```

## 作业标记和元数据

### 添加职位标签

```python
job = dxpy.DXApplet("applet-xxxx").run(
    {"input": "..."},
    tags=["experiment1", "batch2", "production"]
)
```

### 添加作业属性

```python
job = dxpy.DXApplet("applet-xxxx").run(
    {"input": "..."},
    properties={
        "experiment": "exp001",
        "sample": "sample1",
        "batch": "batch2"
    }
)
```

### 寻找工作

```python
# Find jobs by tag
jobs = dxpy.find_jobs(
    project="project-xxxx",
    tags=["experiment1"],
    describe=True
)

for job in jobs:
    print(f"{job['describe']['name']}: {job['id']}")
```

## 最佳实践

1. **作业命名**：使用描述性名称以便于跟踪
2. **标签和属性**：为组织和可搜索性标记作业
3. **资源选择**：为工作负载选择合适的实例类型
4. **错误处理**：检查作业状态并优雅地处理失败
5. **并行处理**：使用子作业进行独立的并行任务
6. **工作流程**：使用工作流程进行复杂的多步骤分析
7. **监控**：监控长时间运行的作业并检查日志是否存在问题
8. **成本管理**：使用适当的实例类型来平衡成本/性能
9. **超时**：设置合理的超时，防止作业失控
10. **清理**：删除失败或过时的作业

## 调试技巧

1. **检查日志**：始终查看作业日志中的错误消息
2. **验证输入**：确保输入文件已关闭且可访问
3. **本地测试**：部署到平台之前在本地测试逻辑
4. **从小规模开始**：在扩大规模之前使用小数据集进行测试
5. **监控资源**：检查作业是否耗尽内存或磁盘空间
6. **实例类型**：如果作业因资源原因失败，请尝试更大的实例