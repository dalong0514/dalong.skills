<!-- 此文件由机器翻译自 api_examples.md -->

# UniProt API 示例

以多种语言与 UniProt REST API 交互的实用代码示例。

## Python 示例

### 示例 1：基本搜索
```python
import requests

# Search for human insulin proteins
url = "https://rest.uniprot.org/uniprotkb/search"
params = {
    "query": "insulin AND organism_id:9606 AND reviewed:true",
    "format": "json",
    "size": 10
}

response = requests.get(url, params=params)
data = response.json()

for result in data['results']:
    print(f"{result['primaryAccession']}: {result['proteinDescription']['recommendedName']['fullName']['value']}")
```

### 示例 2：检索蛋白质序列
<<<代码块_1>>>

### 示例 3：自定义字段
<<<代码块_2>>>

### 示例 4：ID 映射
<<<代码块_3>>>

### 示例 5：流式传输大型结果
<<<代码块_4>>>

### 示例 6：分页
<<<代码块_5>>>

## 卷曲示例

### 示例 1：简单搜索
<<<代码块_6>>>

### 示例 2：获取蛋白质条目
```bash
# Get human insulin in FASTA format
curl "https://rest.uniprot.org/uniprotkb/P01308.fasta"
```

### 示例 3：自定义字段
```bash
# Get specific fields in TSV format
curl "https://rest.uniprot.org/uniprotkb/search?query=gene:BRCA1&format=tsv&fields=accession,gene_names,length"
```

### 示例 4：ID 映射 - 提交作业
```bash
# Submit mapping job
curl -X POST "https://rest.uniprot.org/idmapping/run" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "from=UniProtKB_AC-ID&to=PDB&ids=P01308,P04637"
```

### 示例 5：ID 映射 - 获取结果
```bash
# Get mapping results (replace JOB_ID)
curl "https://rest.uniprot.org/idmapping/results/JOB_ID"
```

### 示例 6：下载所有结果
```bash
# Download all human reviewed proteins
curl "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta" \
  -o human_proteins.fasta
```

## R 示例

### 示例 1：基本搜索
```r
library(httr)
library(jsonlite)

# Search for insulin proteins
url <- "https://rest.uniprot.org/uniprotkb/search"
query_params <- list(
  query = "insulin AND organism_id:9606",
  format = "json",
  size = 10
)

response <- GET(url, query = query_params)
data <- fromJSON(content(response, "text"))

# Extract accessions and names
proteins <- data$results[, c("primaryAccession", "proteinDescription")]
print(proteins)
```

### 示例 2：获取序列
```r
library(httr)

# Get protein sequence
accession <- "P01308"
url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".fasta")

response <- GET(url)
sequence <- content(response, "text")
cat(sequence)
```

### 示例 3：下载到数据框
```r
library(httr)
library(readr)

# Get data as TSV
url <- "https://rest.uniprot.org/uniprotkb/search"
query_params <- list(
  query = "gene:BRCA1 AND reviewed:true",
  format = "tsv",
  fields = "accession,gene_names,organism_name,length"
)

response <- GET(url, query = query_params)
data <- read_tsv(content(response, "text"))
print(data)
```

## JavaScript 示例

### 示例 1：获取 API
```javascript
// Search for proteins
async function searchUniProt(query) {
  const url = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(query)}&format=json&size=10`;

  const response = await fetch(url);
  const data = await response.json();

  return data.results;
}

// Usage
searchUniProt("insulin AND organism_id:9606")
  .then(results => console.log(results));
```

### 示例 2：获取蛋白质条目
```javascript
async function getProtein(accession, format = "json") {
  const url = `https://rest.uniprot.org/uniprotkb/${accession}.${format}`;

  const response = await fetch(url);

  if (format === "json") {
    return await response.json();
  } else {
    return await response.text();
  }
}

// Usage
getProtein("P01308", "fasta")
  .then(sequence => console.log(sequence));
```

### 示例 3：ID 映射
```javascript
async function mapIds(ids, fromDb, toDb) {
  // Submit job
  const submitUrl = "https://rest.uniprot.org/idmapping/run";
  const formData = new URLSearchParams({
    from: fromDb,
    to: toDb,
    ids: ids.join(",")
  });

  const submitResponse = await fetch(submitUrl, {
    method: "POST",
    body: formData
  });
  const { jobId } = await submitResponse.json();

  // Poll for completion
  const statusUrl = `https://rest.uniprot.org/idmapping/status/${jobId}`;
  while (true) {
    const statusResponse = await fetch(statusUrl);
    const status = await statusResponse.json();

    if ("results" in status || "failedIds" in status) {
      break;
    }

    await new Promise(resolve => setTimeout(resolve, 3000));
  }

  // Get results
  const resultsUrl = `https://rest.uniprot.org/idmapping/results/${jobId}`;
  const resultsResponse = await fetch(resultsUrl);
  return await resultsResponse.json();
}

// Usage
mapIds(["P01308", "P04637"], "UniProtKB_AC-ID", "PDB")
  .then(mapping => console.log(mapping));
```

## 高级示例

### 示例：具有速率限制的批处理
```python
import requests
import time
from typing import List, Dict

class UniProtClient:
    def __init__(self, rate_limit=1.0):
        self.base_url = "https://rest.uniprot.org"
        self.rate_limit = rate_limit
        self.last_request = 0

    def _rate_limit(self):
        """Enforce rate limiting"""
        elapsed = time.time() - self.last_request
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        self.last_request = time.time()

    def batch_get_proteins(self, accessions: List[str],
                          batch_size: int = 100) -> List[Dict]:
        """Get proteins in batches"""
        results = []

        for i in range(0, len(accessions), batch_size):
            batch = accessions[i:i + batch_size]
            query = " OR ".join([f"accession:{acc}" for acc in batch])

            self._rate_limit()

            response = requests.get(
                f"{self.base_url}/uniprotkb/search",
                params={
                    "query": query,
                    "format": "json",
                    "size": batch_size
                }
            )

            if response.ok:
                data = response.json()
                results.extend(data.get('results', []))
            else:
                print(f"Error in batch {i//batch_size}: {response.status_code}")

        return results

# Usage
client = UniProtClient(rate_limit=0.5)
accessions = ["P01308", "P04637", "P12345", "Q9Y6K9"]
proteins = client.batch_get_proteins(accessions)
```

### 示例：带进度条下载
```python
import requests
from tqdm import tqdm

def download_with_progress(query, output_file, format="fasta"):
    """Download results with progress bar"""
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": query,
        "format": format
    }

    response = requests.get(url, params=params, stream=True)
    total_size = int(response.headers.get('content-length', 0))

    with open(output_file, 'wb') as f, \
         tqdm(total=total_size, unit='B', unit_scale=True) as pbar:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
            pbar.update(len(chunk))

# Usage
download_with_progress(
    "organism_id:9606 AND reviewed:true",
    "human_proteome.fasta"
)
```

## 资源

- API 文档：https://www.uniprot.org/help/api
- 交互式 API 资源管理器：https://www.uniprot.org/api-documentation
- Python 客户端（未压缩）：https://github.com/multimeric/Unipressed
- 生物服务包：https://bioservices.readthedocs.io/