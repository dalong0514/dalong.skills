<!-- 此文件由机器翻译自 scheduled-jobs.md -->

# 预定的作业和 Cron

## 基本调度

安排功能定期或特定时间自动运行。

### 简单的每日时间表

```python
import modal

app = modal.App()

@app.function(schedule=modal.Period(days=1))
def daily_task():
    print("Running daily task")
    # Process data, send reports, etc.
```

部署激活：
<<<代码块_1>>>

函数从部署时间起每 24 小时运行一次。

## 时间表类型

### 期间时间表

从部署时间开始以固定间隔运行：

<<<代码块_2>>>

**注意**：重新部署会重置周期计时器。

### Cron 计划

使用 cron 语法在特定时间运行：

<<<代码块_3>>>

**Cron 语法**：`minute hour day month day_of_week`
- 分钟：0-59
- 时间：0-23
- 日：1-31
- 月份：1-12
- 星期几：0-6（0 = 星期日）

### 时区支持

为 cron 计划指定时区：

<<<代码块_4>>>

## 部署

### 部署预定功能

<<<代码块_5>>>

计划的功能将持续存在，直到明确停止为止。

### 程序化部署

<<<代码块_6>>>

## 监控

### 查看执行日志

检查 https://modal.com/apps 是否有：
- 过去的执行日志
- 执行历史
- 失败通知

### 手动运行

通过仪表板“立即运行”按钮立即触发预定功能。

## 日程管理

### 暂停日程

计划无法暂停。停止：
1. 删除`schedule`参数
2. 重新部署应用程序

### 更新时间表

更改计划参数并重新部署：

```python
# Update from daily to weekly
@app.function(schedule=modal.Period(days=7))
def task():
    ...
```

```bash
modal deploy script.py
```

## 常见模式

### 数据管道

```python
@app.function(
    schedule=modal.Cron("0 2 * * *"),  # 2 AM daily
    timeout=3600,                       # 1 hour timeout
)
def etl_pipeline():
    # Extract data from sources
    data = extract_data()

    # Transform data
    transformed = transform_data(data)

    # Load to warehouse
    load_to_warehouse(transformed)
```

### 模型再训练

```python
volume = modal.Volume.from_name("models")

@app.function(
    schedule=modal.Cron("0 0 * * 0"),  # Weekly on Sunday midnight
    gpu="A100",
    timeout=7200,                       # 2 hours
    volumes={"/models": volume}
)
def retrain_model():
    # Load latest data
    data = load_training_data()

    # Train model
    model = train(data)

    # Save new model
    save_model(model, "/models/latest.pt")
    volume.commit()
```

### 报告生成

```python
@app.function(
    schedule=modal.Cron("0 9 * * 1"),  # Monday 9 AM
    secrets=[modal.Secret.from_name("email-creds")]
)
def weekly_report():
    # Generate report
    report = generate_analytics_report()

    # Send email
    send_email(
        to="team@company.com",
        subject="Weekly Analytics Report",
        body=report
    )
```

### 数据清理

```python
@app.function(schedule=modal.Period(hours=6))
def cleanup_old_data():
    # Remove data older than 30 days
    cutoff = datetime.now() - timedelta(days=30)
    delete_old_records(cutoff)
```

## 配置机密和卷

预定函数支持所有函数参数：

```python
vol = modal.Volume.from_name("data")
secret = modal.Secret.from_name("api-keys")

@app.function(
    schedule=modal.Cron("0 */6 * * *"),  # Every 6 hours
    secrets=[secret],
    volumes={"/data": vol},
    cpu=4.0,
    memory=16384,
)
def sync_data():
    import os

    api_key = os.environ["API_KEY"]

    # Fetch from external API
    data = fetch_external_data(api_key)

    # Save to volume
    with open("/data/latest.json", "w") as f:
        json.dump(data, f)

    vol.commit()
```

## 动态调度

以编程方式更新时间表：

```python
@app.function()
def main_task():
    ...

@app.function(schedule=modal.Cron("0 6 * * *", timezone="America/New_York"))
def enable_high_traffic_mode():
    main_task.update_autoscaler(min_containers=5)

@app.function(schedule=modal.Cron("0 22 * * *", timezone="America/New_York"))
def disable_high_traffic_mode():
    main_task.update_autoscaler(min_containers=0)
```

## 错误处理

失败的预定功能将：
- 在仪表板中显示故障
- 发送通知（可配置）
- 在下一次计划运行时重试

```python
@app.function(
    schedule=modal.Cron("0 * * * *"),
    retries=3,  # Retry failed runs
    timeout=1800
)
def robust_task():
    try:
        perform_task()
    except Exception as e:
        # Log error
        print(f"Task failed: {e}")
        # Optionally send alert
        send_alert(f"Scheduled task failed: {e}")
        raise
```

## 最佳实践

1. **设置超时**：始终为预定功能指定超时
2. **使用适当的时间表**：Period 表示相对时间，Cron 表示绝对时间
3. **监控故障**：定期检查仪表板是否有失败的运行
4. **幂等操作**：设计任务以安全地处理重新运行
5. **资源限制**：为计划的工作负载设置适当的CPU/内存
6. **时区感知**：为 cron 计划指定时区