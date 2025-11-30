<!-- 此文件由机器翻译自 examples.md -->

# 科学计算的常见模式

## 机器学习模型推理

### 基本模型服务

```python
import modal

app = modal.App("ml-inference")

image = (
    modal.Image.debian_slim()
    .uv_pip_install("torch", "transformers")
)

@app.cls(
    image=image,
    gpu="L40S",
)
class Model:
    @modal.enter()
    def load_model(self):
        from transformers import AutoModel, AutoTokenizer
        self.model = AutoModel.from_pretrained("bert-base-uncased")
        self.tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

    @modal.method()
    def predict(self, text: str):
        inputs = self.tokenizer(text, return_tensors="pt")
        outputs = self.model(**inputs)
        return outputs.last_hidden_state.mean(dim=1).tolist()

@app.local_entrypoint()
def main():
    model = Model()
    result = model.predict.remote("Hello world")
    print(result)
```

### 批量服务模型

<<<代码块_1>>>

## 批处理

### 并行数据处理

<<<代码块_2>>>

### 有进度的批处理

<<<代码块_3>>>

## 数据分析管道

### ETL 管道

<<<代码块_4>>>

## GPU 加速计算

### 分布式训练

<<<代码块_5>>>

### GPU 批量推理

<<<代码块_6>>>

## 科学计算

### 分子动力学模拟

```python
@app.function(
    image=modal.Image.debian_slim().apt_install("openmpi-bin").uv_pip_install("mpi4py", "numpy"),
    cpu=16.0,
    memory=65536,
    timeout=7200,
)
def run_simulation(config: dict):
    import numpy as np

    # Initialize system
    positions = initialize_positions(config["n_particles"])
    velocities = initialize_velocities(config["temperature"])

    # Run MD steps
    for step in range(config["n_steps"]):
        forces = compute_forces(positions)
        velocities += forces * config["dt"]
        positions += velocities * config["dt"]

        if step % 1000 == 0:
            energy = compute_energy(positions, velocities)
            print(f"Step {step}, Energy: {energy}")

    return positions, velocities
```

### 分布式蒙特卡罗

```python
@app.function(cpu=2.0)
def monte_carlo_trial(trial_id: int, n_samples: int):
    import random

    count = sum(1 for _ in range(n_samples)
                if random.random()**2 + random.random()**2 <= 1)

    return count

@app.local_entrypoint()
def estimate_pi():
    n_trials = 100
    n_samples_per_trial = 1_000_000

    # Run trials in parallel
    results = list(monte_carlo_trial.map(
        range(n_trials),
        [n_samples_per_trial] * n_trials
    ))

    total_count = sum(results)
    total_samples = n_trials * n_samples_per_trial

    pi_estimate = 4 * total_count / total_samples
    print(f"Estimated π = {pi_estimate}")
```

## 卷数据处理

### 图像处理管道

```python
volume = modal.Volume.from_name("images")
IMAGE_PATH = "/images"

@app.function(
    image=modal.Image.debian_slim().uv_pip_install("Pillow", "numpy"),
    volumes={IMAGE_PATH: volume}
)
def process_image(filename: str):
    from PIL import Image
    import numpy as np

    # Load image
    img = Image.open(f"{IMAGE_PATH}/raw/{filename}")

    # Process
    img_array = np.array(img)
    processed = apply_filters(img_array)

    # Save
    result_img = Image.fromarray(processed)
    result_img.save(f"{IMAGE_PATH}/processed/{filename}")

    return filename

@app.function(volumes={IMAGE_PATH: volume})
def process_all_images():
    import os

    # Get all images
    filenames = os.listdir(f"{IMAGE_PATH}/raw")

    # Process in parallel
    results = list(process_image.map(filenames))

    volume.commit()
    return f"Processed {len(results)} images"
```

## 用于科学计算的 Web API

```python
image = modal.Image.debian_slim().uv_pip_install("fastapi[standard]", "numpy", "scipy")

@app.function(image=image)
@modal.fastapi_endpoint(method="POST")
def compute_statistics(data: dict):
    import numpy as np
    from scipy import stats

    values = np.array(data["values"])

    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values)),
        "skewness": float(stats.skew(values)),
        "kurtosis": float(stats.kurtosis(values))
    }
```

## 预定数据收集

```python
@app.function(
    schedule=modal.Cron("*/30 * * * *"),  # Every 30 minutes
    secrets=[modal.Secret.from_name("api-keys")],
    volumes={"/data": modal.Volume.from_name("sensor-data")}
)
def collect_sensor_data():
    import requests
    import json
    from datetime import datetime

    # Fetch from API
    response = requests.get(
        "https://api.example.com/sensors",
        headers={"Authorization": f"Bearer {os.environ['API_KEY']}"}
    )

    data = response.json()

    # Save with timestamp
    timestamp = datetime.now().isoformat()
    with open(f"/data/{timestamp}.json", "w") as f:
        json.dump(data, f)

    volume.commit()

    return f"Collected {len(data)} sensor readings"
```

## 最佳实践

### 使用类来处理有状态的工作负载

```python
@app.cls(gpu="A100")
class ModelService:
    @modal.enter()
    def setup(self):
        # Load once, reuse across requests
        self.model = load_heavy_model()

    @modal.method()
    def predict(self, x):
        return self.model(x)
```

### 批量相似的工作负载

```python
@app.function()
def process_many(items: list):
    # More efficient than processing one at a time
    return [process(item) for item in items]
```

### 对大型数据集使用卷

```python
# Store large datasets in volumes, not in image
volume = modal.Volume.from_name("dataset")

@app.function(volumes={"/data": volume})
def train():
    data = load_from_volume("/data/training.parquet")
    model = train_model(data)
```

### 扩展到 GPU 之前的配置文件

```python
# Test on CPU first
@app.function(cpu=4.0)
def test_pipeline():
    ...

# Then scale to GPU if needed
@app.function(gpu="A100")
def gpu_pipeline():
    ...
```