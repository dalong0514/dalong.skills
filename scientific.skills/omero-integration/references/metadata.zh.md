<!-- 此文件由机器翻译自 metadata.md -->

# 元数据和注释

本参考涵盖在 OMERO 中创建和管理注释，包括标签、键值对、文件附件和注释。

## 注释类型

OMERO 支持多种注释类型：

- **TagAnnotation**：用于分类的文本标签
- **MapAnnotation**：结构化元数据的键值对
- **FileAnnotation**：文件附件（PDF、CSV、分析结果等）
- **CommentAnnotation**：自由文本评论
- **LongAnnotation**：整数值
- **DoubleAnnotation**：浮点值
- **BooleanAnnotation**：布尔值
- **TimestampAnnotation**：日期/时间戳
- **术语注释**：本体术语

## 标签注释

### 创建并链接标签

```python
import omero.gateway

# Create new tag
tag_ann = omero.gateway.TagAnnotationWrapper(conn)
tag_ann.setValue("Experiment 2024")
tag_ann.setDescription("Optional description of this tag")
tag_ann.save()

# Link tag to an object
project = conn.getObject("Project", project_id)
project.linkAnnotation(tag_ann)
```

### 使用命名空间创建标签

<<<代码块_1>>>

### 重用现有标签

<<<代码块_2>>>

### 创建标签集（带子标签）

<<<代码块_3>>>

## 映射注释（键值对）

### 创建地图注释

<<<代码块_4>>>

### 地图注释的自定义命名空间

<<<代码块_5>>>

### 读取地图注释

<<<代码块_6>>>

## 文件注释

### 上传并附加文件

```python
import os

# Prepare file
file_path = "analysis_results.csv"

# Create file annotation
namespace = "mylab.analysis.results"
file_ann = conn.createFileAnnfromLocalFile(
    file_path,
    mimetype="text/csv",
    ns=namespace,
    desc="Cell segmentation results"
)

# Link to dataset
dataset = conn.getObject("Dataset", dataset_id)
dataset.linkAnnotation(file_ann)
```

### 支持的 MIME 类型

常见的 MIME 类型：
- 文本：`"text/plain"`、`"text/csv"`、`"text/tab-separated-values"`
- 文档：`"application/pdf"`、`"application/vnd.ms-excel"`
- 图片：`"image/png"`、`"image/jpeg"`
- 数据：`"application/json"`、`"application/xml"`
- 档案：`"application/zip"`、`"application/gzip"`

### 上传多个文件

```python
files = ["figure1.pdf", "figure2.pdf", "table1.csv"]
namespace = "publication.supplementary"

dataset = conn.getObject("Dataset", dataset_id)

for file_path in files:
    file_ann = conn.createFileAnnfromLocalFile(
        file_path,
        mimetype="application/octet-stream",
        ns=namespace,
        desc=f"Supplementary file: {os.path.basename(file_path)}"
    )
    dataset.linkAnnotation(file_ann)
```

### 下载文件注释

```python
import os

# Get object with file annotation
image = conn.getObject("Image", image_id)

# Download directory
download_path = "./downloads"
os.makedirs(download_path, exist_ok=True)

# Filter by namespace
namespace = "mylab.analysis.results"

for ann in image.listAnnotations(ns=namespace):
    if isinstance(ann, omero.gateway.FileAnnotationWrapper):
        file_name = ann.getFile().getName()
        file_path = os.path.join(download_path, file_name)

        print(f"Downloading: {file_name}")

        # Download file in chunks
        with open(file_path, 'wb') as f:
            for chunk in ann.getFileInChunks():
                f.write(chunk)

        print(f"Saved to: {file_path}")
```

### 获取文件注释元数据

```python
for ann in dataset.listAnnotations():
    if isinstance(ann, omero.gateway.FileAnnotationWrapper):
        orig_file = ann.getFile()

        print(f"File Annotation ID: {ann.getId()}")
        print(f"  File Name: {orig_file.getName()}")
        print(f"  File Size: {orig_file.getSize()} bytes")
        print(f"  MIME Type: {orig_file.getMimetype()}")
        print(f"  Namespace: {ann.getNs()}")
        print(f"  Description: {ann.getDescription()}")
```

## 注释注释

### 添加评论

```python
# Create comment
comment = omero.gateway.CommentAnnotationWrapper(conn)
comment.setValue("This image shows excellent staining quality")
comment.save()

# Link to image
image = conn.getObject("Image", image_id)
image.linkAnnotation(comment)
```

### 添加带有命名空间的注释

```python
comment = omero.gateway.CommentAnnotationWrapper(conn)
comment.setValue("Approved for publication")
comment.setNs("mylab.publication.status")
comment.save()

dataset = conn.getObject("Dataset", dataset_id)
dataset.linkAnnotation(comment)
```

## 数字注释

### 长注释（整数）

```python
# Create long annotation
long_ann = omero.gateway.LongAnnotationWrapper(conn)
long_ann.setValue(42)
long_ann.setNs("mylab.cell.count")
long_ann.save()

image = conn.getObject("Image", image_id)
image.linkAnnotation(long_ann)
```

### 双注释（浮动）

```python
# Create double annotation
double_ann = omero.gateway.DoubleAnnotationWrapper(conn)
double_ann.setValue(3.14159)
double_ann.setNs("mylab.fluorescence.intensity")
double_ann.save()

image = conn.getObject("Image", image_id)
image.linkAnnotation(double_ann)
```

## 列出注释

### 列出对象上的所有注释

```python
import omero.model

# Get object
project = conn.getObject("Project", project_id)

# List all annotations
for ann in project.listAnnotations():
    print(f"Annotation ID: {ann.getId()}")
    print(f"  Type: {ann.OMERO_TYPE}")
    print(f"  Added by: {ann.link.getDetails().getOwner().getOmeName()}")

    # Type-specific handling
    if ann.OMERO_TYPE == omero.model.TagAnnotationI:
        print(f"  Tag value: {ann.getTextValue()}")

    elif isinstance(ann, omero.gateway.MapAnnotationWrapper):
        print(f"  Map data: {ann.getValue()}")

    elif isinstance(ann, omero.gateway.FileAnnotationWrapper):
        print(f"  File: {ann.getFile().getName()}")

    elif isinstance(ann, omero.gateway.CommentAnnotationWrapper):
        print(f"  Comment: {ann.getValue()}")

    print()
```

### 按命名空间过滤注释

```python
# Get annotations with specific namespace
namespace = "mylab.qc.tags"

for ann in image.listAnnotations(ns=namespace):
    print(f"Found annotation: {ann.getId()}")

    if isinstance(ann, omero.gateway.MapAnnotationWrapper):
        for key, value in ann.getValue():
            print(f"  {key}: {value}")
```

### 获取第一个带有命名空间的注释

```python
# Get single annotation by namespace
namespace = "mylab.analysis.results"
ann = dataset.getAnnotation(namespace)

if ann:
    print(f"Found annotation with namespace: {ann.getNs()}")
else:
    print("No annotation found with that namespace")
```

### 跨多个对象查询注释

```python
# Get all tag annotations linked to image IDs
image_ids = [1, 2, 3, 4, 5]

for link in conn.getAnnotationLinks('Image', parent_ids=image_ids):
    ann = link.getChild()

    if isinstance(ann._obj, omero.model.TagAnnotationI):
        print(f"Image {link.getParent().getId()}: Tag '{ann.getTextValue()}'")
```

## 计算注释

```python
# Count annotations on project
project_id = 123
count = conn.countAnnotations('Project', [project_id])
print(f"Project has {count[project_id]} annotations")

# Count annotations on multiple images
image_ids = [1, 2, 3]
counts = conn.countAnnotations('Image', image_ids)

for image_id, count in counts.items():
    print(f"Image {image_id}: {count} annotations")
```

## 注释链接

### 手动创建注释链接

```python
# Get annotation and image
tag = conn.getObject("TagAnnotation", tag_id)
image = conn.getObject("Image", image_id)

# Create link
link = omero.model.ImageAnnotationLinkI()
link.setParent(omero.model.ImageI(image.getId(), False))
link.setChild(omero.model.TagAnnotationI(tag.getId(), False))

# Save link
conn.getUpdateService().saveAndReturnObject(link)
```

### 更新注释链接

```python
# Get existing links
annotation_ids = [1, 2, 3]
new_tag_id = 5

for link in conn.getAnnotationLinks('Image', ann_ids=annotation_ids):
    print(f"Image ID: {link.getParent().id}")

    # Change linked annotation
    link._obj.child = omero.model.TagAnnotationI(new_tag_id, False)
    link.save()
```

## 删除注释

### 删除注释

```python
# Get image
image = conn.getObject("Image", image_id)

# Collect annotation IDs to delete
to_delete = []
namespace = "mylab.temp.annotations"

for ann in image.listAnnotations(ns=namespace):
    to_delete.append(ann.getId())

# Delete annotations
if to_delete:
    conn.deleteObjects('Annotation', to_delete, wait=True)
    print(f"Deleted {len(to_delete)} annotations")
```

### 取消链接注释（保留注释、删除链接）

```python
# Get image
image = conn.getObject("Image", image_id)

# Collect link IDs to delete
to_delete = []

for ann in image.listAnnotations():
    if isinstance(ann, omero.gateway.TagAnnotationWrapper):
        to_delete.append(ann.link.getId())

# Delete links (annotations remain in database)
if to_delete:
    conn.deleteObjects("ImageAnnotationLink", to_delete, wait=True)
    print(f"Unlinked {len(to_delete)} annotations")
```

### 删除特定注释类型

```python
import omero.gateway

# Delete only map annotations
image = conn.getObject("Image", image_id)
to_delete = []

for ann in image.listAnnotations():
    if isinstance(ann, omero.gateway.MapAnnotationWrapper):
        to_delete.append(ann.getId())

conn.deleteObjects('Annotation', to_delete, wait=True)
```

## 注释所有权

### 设置注释所有者（仅限管理员）

```python
import omero.model

# Create tag with specific owner
tag_ann = omero.gateway.TagAnnotationWrapper(conn)
tag_ann.setValue("Admin Tag")

# Set owner (requires admin privileges)
user_id = 5
tag_ann._obj.details.owner = omero.model.ExperimenterI(user_id, False)
tag_ann.save()
```

### 以另一个用户身份创建注释（仅限管理员）

```python
# Admin connection
admin_conn = BlitzGateway(admin_user, admin_pass, host=host, port=4064)
admin_conn.connect()

# Get target user
user_id = 10
user = admin_conn.getObject("Experimenter", user_id).getName()

# Create connection as user
user_conn = admin_conn.suConn(user)

# Create annotation as that user
map_ann = omero.gateway.MapAnnotationWrapper(user_conn)
map_ann.setNs("mylab.metadata")
map_ann.setValue([["key", "value"]])
map_ann.save()

# Link to project
project = admin_conn.getObject("Project", project_id)
project.linkAnnotation(map_ann)

# Close connections
user_conn.close()
admin_conn.close()
```

## 批量标注操作

### 标记多个图像

```python
# Create or get tag
tag = omero.gateway.TagAnnotationWrapper(conn)
tag.setValue("Validated")
tag.save()

# Get images to tag
dataset = conn.getObject("Dataset", dataset_id)

# Tag all images in dataset
for image in dataset.listChildren():
    image.linkAnnotation(tag)
    print(f"Tagged image: {image.getName()}")
```

### 批量添加地图标注

```python
# Prepare metadata for multiple images
image_metadata = {
    101: [["Quality", "Good"], ["Reviewed", "Yes"]],
    102: [["Quality", "Excellent"], ["Reviewed", "Yes"]],
    103: [["Quality", "Poor"], ["Reviewed", "No"]]
}

# Add annotations
for image_id, kv_data in image_metadata.items():
    image = conn.getObject("Image", image_id)

    if image:
        map_ann = omero.gateway.MapAnnotationWrapper(conn)
        map_ann.setNs("mylab.qc")
        map_ann.setValue(kv_data)
        map_ann.save()

        image.linkAnnotation(map_ann)
        print(f"Annotated image {image_id}")
```

## 命名空间

### 标准 OMERO 命名空间

```python
import omero.constants.metadata as omero_ns

# Client map annotation namespace
omero_ns.NSCLIENTMAPANNOTATION
# "openmicroscopy.org/omero/client/mapAnnotation"

# Bulk annotations namespace
omero_ns.NSBULKANNOTATIONS
# "openmicroscopy.org/omero/bulk_annotations"
```

### 自定义命名空间

自定义命名空间的最佳实践：
- 使用反向域表示法：`"org.mylab.category.subcategory"`
- 具体：`"com.company.project.analysis.v1"`
- 如果架构可能更改，请包含版本：`"mylab.metadata.v2"`

```python
# Define namespaces
NS_QC = "org.mylab.quality_control"
NS_ANALYSIS = "org.mylab.image_analysis.v1"
NS_PUBLICATION = "org.mylab.publication.2024"

# Use in annotations
map_ann.setNs(NS_ANALYSIS)
```

## 按类型加载所有注释

### 加载所有文件注释

```python
# Define namespaces to include/exclude
ns_to_include = ["mylab.analysis.results"]
ns_to_exclude = []

# Get metadata service
metadataService = conn.getMetadataService()

# Load all file annotations with namespace
annotations = metadataService.loadSpecifiedAnnotations(
    'omero.model.FileAnnotation',
    ns_to_include,
    ns_to_exclude,
    None
)

for ann in annotations:
    print(f"File Annotation ID: {ann.getId().getValue()}")
    print(f"  File: {ann.getFile().getName().getValue()}")
    print(f"  Size: {ann.getFile().getSize().getValue()} bytes")
```

## 完整示例

```python
from omero.gateway import BlitzGateway
import omero.gateway
import omero.constants.metadata

HOST = 'omero.example.com'
PORT = 4064
USERNAME = 'user'
PASSWORD = 'pass'

with BlitzGateway(USERNAME, PASSWORD, host=HOST, port=PORT) as conn:
    # Get dataset
    dataset = conn.getObject("Dataset", dataset_id)

    # Add tag
    tag = omero.gateway.TagAnnotationWrapper(conn)
    tag.setValue("Analysis Complete")
    tag.save()
    dataset.linkAnnotation(tag)

    # Add map annotation with metadata
    metadata = [
        ["Analysis Date", "2024-10-20"],
        ["Software", "CellProfiler 4.2"],
        ["Pipeline", "cell_segmentation_v3"]
    ]
    map_ann = omero.gateway.MapAnnotationWrapper(conn)
    map_ann.setNs(omero.constants.metadata.NSCLIENTMAPANNOTATION)
    map_ann.setValue(metadata)
    map_ann.save()
    dataset.linkAnnotation(map_ann)

    # Add file annotation
    file_ann = conn.createFileAnnfromLocalFile(
        "analysis_summary.pdf",
        mimetype="application/pdf",
        ns="mylab.reports",
        desc="Analysis summary report"
    )
    dataset.linkAnnotation(file_ann)

    # Add comment
    comment = omero.gateway.CommentAnnotationWrapper(conn)
    comment.setValue("Dataset ready for review")
    comment.save()
    dataset.linkAnnotation(comment)

    print(f"Added 4 annotations to dataset {dataset.getName()}")
```

## 最佳实践

1. **使用命名空间**：始终使用命名空间来组织注释
2. **描述性标签**：使用清晰、一致的标签名称
3. **结构化元数据**：对于结构化数据，更喜欢地图注释而不是注释
4. **文件组织**：使用描述性文件名和 MIME 类型
5. **链接重用**：重用现有标签而不是创建重复项
6. **批量操作**：循环处理多个对象以提高效率
7. **错误处理**：在链接之前检查是否成功保存
8. **Cleanup**：不再需要时删除临时注释
9. **文档**：记录自定义命名空间的含义
10. **权限**：考虑协作工作流程的注释所有权