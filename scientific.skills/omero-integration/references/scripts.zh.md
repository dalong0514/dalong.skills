<!-- 此文件由机器翻译自 scripts.md -->

# 脚本和批量操作

本参考内容涵盖了为服务器端处理和批处理操作创建 OMERO.scripts。

## OMERO.scripts 概述

OMERO.scripts 是在 OMERO 服务器上运行的 Python 脚本，可以从 OMERO 客户端（Web、insight、CLI）调用。它们充当扩展 OMERO 功能的插件。

### 主要特点

- **服务器端执行**：脚本在服务器上运行，避免数据传输
- **客户端集成**：可通过自动生成的 UI 从任何 OMERO 客户端调用
- **参数处理**：通过验证定义输入参数
- **结果报告**：向客户返回图像、文件或消息
- **批处理**：高效处理多个图像或数据集

## 基本脚本结构

### 最小脚本模板

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import omero
from omero.gateway import BlitzGateway
import omero.scripts as scripts
from omero.rtypes import rlong, rstring, robject

def run_script():
    """
    Main script function.
    """
    # Script definition
    client = scripts.client(
        'Script_Name.py',
        """
        Description of what this script does.
        """,

        # Input parameters
        scripts.String("Data_Type", optional=False, grouping="1",
                      description="Choose source of images",
                      values=[rstring('Dataset'), rstring('Image')],
                      default=rstring('Dataset')),

        scripts.Long("IDs", optional=False, grouping="2",
                    description="Dataset or Image ID(s)").ofType(rlong(0)),

        # Outputs
        namespaces=[omero.constants.namespaces.NSDYNAMIC],
        version="1.0"
    )

    try:
        # Get connection
        conn = BlitzGateway(client_obj=client)

        # Get script parameters
        script_params = client.getInputs(unwrap=True)
        data_type = script_params["Data_Type"]
        ids = script_params["IDs"]

        # Process data
        message = process_data(conn, data_type, ids)

        # Return results
        client.setOutput("Message", rstring(message))

    finally:
        client.closeSession()

def process_data(conn, data_type, ids):
    """
    Process images based on parameters.
    """
    # Implementation here
    return "Processing complete"

if __name__ == "__main__":
    run_script()
```

## 脚本参数

### 参数类型

<<<代码块_1>>>

### 参数分组

<<<代码块_2>>>

## 访问输入数据

### 获取脚本参数

<<<代码块_3>>>

### 从参数获取图像

<<<代码块_4>>>

## 处理图像

### 批量图像处理

<<<代码块_5>>>

## 生成输出

### 返回消息

<<<代码块_6>>>

### 返回图像

```python
# Return newly created image
new_image = conn.createImageFromNumpySeq(...)
client.setOutput("New_Image", robject(new_image._obj))
```

### 返回文件

```python
# Create and return file annotation
file_ann = conn.createFileAnnfromLocalFile(
    output_file_path,
    mimetype="text/csv",
    ns="analysis.results"
)

client.setOutput("Result_File", robject(file_ann._obj))
```

### 返回表

```python
# Create OMERO table and return
resources = conn.c.sf.sharedResources()
table = create_results_table(resources, results)
orig_file = table.getOriginalFile()
table.close()

# Create file annotation
file_ann = omero.model.FileAnnotationI()
file_ann.setFile(orig_file)
file_ann = conn.getUpdateService().saveAndReturnObject(file_ann)

client.setOutput("Results_Table", robject(file_ann._obj))
```

## 完整的示例脚本

### 示例 1：最大强度投影

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import omero
from omero.gateway import BlitzGateway
import omero.scripts as scripts
from omero.rtypes import rlong, rstring, robject
import numpy as np

def run_script():
    client = scripts.client(
        'Maximum_Intensity_Projection.py',
        """
        Creates maximum intensity projection from Z-stack images.
        """,

        scripts.String("Data_Type", optional=False, grouping="1",
                      description="Process images from",
                      values=[rstring('Dataset'), rstring('Image')],
                      default=rstring('Image')),

        scripts.List("IDs", optional=False, grouping="2",
                    description="Dataset or Image ID(s)").ofType(rlong(0)),

        scripts.Bool("Link_to_Source", optional=True, grouping="3",
                    description="Link results to source dataset",
                    default=True),

        version="1.0"
    )

    try:
        conn = BlitzGateway(client_obj=client)
        script_params = client.getInputs(unwrap=True)

        # Get images
        images = get_images(conn, script_params)
        created_images = []

        for image in images:
            print(f"Processing: {image.getName()}")

            # Create MIP
            mip_image = create_mip(conn, image)
            if mip_image:
                created_images.append(mip_image)

        # Report results
        if created_images:
            message = f"Created {len(created_images)} MIP images"
            # Return first image for display
            client.setOutput("Message", rstring(message))
            client.setOutput("Result", robject(created_images[0]._obj))
        else:
            client.setOutput("Message", rstring("No images created"))

    finally:
        client.closeSession()

def get_images(conn, script_params):
    """Get images from script parameters."""
    images = []
    data_type = script_params["Data_Type"]
    ids = script_params["IDs"]

    if data_type == "Dataset":
        for dataset_id in ids:
            dataset = conn.getObject("Dataset", dataset_id)
            if dataset:
                images.extend(list(dataset.listChildren()))
    else:
        for image_id in ids:
            image = conn.getObject("Image", image_id)
            if image:
                images.append(image)

    return images

def create_mip(conn, source_image):
    """Create maximum intensity projection."""
    pixels = source_image.getPrimaryPixels()
    size_z = source_image.getSizeZ()
    size_c = source_image.getSizeC()
    size_t = source_image.getSizeT()

    if size_z == 1:
        print("  Skipping (single Z-section)")
        return None

    def plane_gen():
        for c in range(size_c):
            for t in range(size_t):
                # Get Z-stack
                z_stack = []
                for z in range(size_z):
                    plane = pixels.getPlane(z, c, t)
                    z_stack.append(plane)

                # Maximum projection
                max_proj = np.max(z_stack, axis=0)
                yield max_proj

    # Create new image
    new_image = conn.createImageFromNumpySeq(
        plane_gen(),
        f"{source_image.getName()}_MIP",
        1, size_c, size_t,
        description="Maximum intensity projection",
        dataset=source_image.getParent()
    )

    return new_image

if __name__ == "__main__":
    run_script()
```

### 示例 2：批量 ROI 分析

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import omero
from omero.gateway import BlitzGateway
import omero.scripts as scripts
from omero.rtypes import rlong, rstring, robject
import omero.grid

def run_script():
    client = scripts.client(
        'Batch_ROI_Analysis.py',
        """
        Analyzes ROIs across multiple images and creates results table.
        """,

        scripts.Long("Dataset_ID", optional=False,
                    description="Dataset with images and ROIs").ofType(rlong(0)),

        scripts.Long("Channel_Index", optional=True,
                    description="Channel to analyze (0-indexed)",
                    default=0, min=0),

        version="1.0"
    )

    try:
        conn = BlitzGateway(client_obj=client)
        script_params = client.getInputs(unwrap=True)

        dataset_id = script_params["Dataset_ID"]
        channel_index = script_params["Channel_Index"]

        # Get dataset
        dataset = conn.getObject("Dataset", dataset_id)
        if not dataset:
            client.setOutput("Message", rstring("Dataset not found"))
            return

        # Analyze ROIs
        results = analyze_rois(conn, dataset, channel_index)

        # Create table
        table_file = create_results_table(conn, dataset, results)

        # Report
        message = f"Analyzed {len(results)} ROIs from {dataset.getName()}"
        client.setOutput("Message", rstring(message))
        client.setOutput("Results_Table", robject(table_file._obj))

    finally:
        client.closeSession()

def analyze_rois(conn, dataset, channel_index):
    """Analyze all ROIs in dataset images."""
    roi_service = conn.getRoiService()
    results = []

    for image in dataset.listChildren():
        result = roi_service.findByImage(image.getId(), None)

        if not result.rois:
            continue

        # Get shape IDs
        shape_ids = []
        for roi in result.rois:
            for shape in roi.copyShapes():
                shape_ids.append(shape.id.val)

        # Get statistics
        stats = roi_service.getShapeStatsRestricted(
            shape_ids, 0, 0, [channel_index]
        )

        # Store results
        for i, stat in enumerate(stats):
            results.append({
                'image_id': image.getId(),
                'image_name': image.getName(),
                'shape_id': shape_ids[i],
                'mean': stat.mean[channel_index],
                'min': stat.min[channel_index],
                'max': stat.max[channel_index],
                'sum': stat.sum[channel_index],
                'area': stat.pointsCount[channel_index]
            })

    return results

def create_results_table(conn, dataset, results):
    """Create OMERO table from results."""
    # Prepare data
    image_ids = [r['image_id'] for r in results]
    shape_ids = [r['shape_id'] for r in results]
    means = [r['mean'] for r in results]
    mins = [r['min'] for r in results]
    maxs = [r['max'] for r in results]
    sums = [r['sum'] for r in results]
    areas = [r['area'] for r in results]

    # Create table
    resources = conn.c.sf.sharedResources()
    repository_id = resources.repositories().descriptions[0].getId().getValue()
    table = resources.newTable(repository_id, f"ROI_Analysis_{dataset.getId()}")

    # Define columns
    columns = [
        omero.grid.ImageColumn('Image', 'Source image', []),
        omero.grid.LongColumn('ShapeID', 'ROI shape ID', []),
        omero.grid.DoubleColumn('Mean', 'Mean intensity', []),
        omero.grid.DoubleColumn('Min', 'Min intensity', []),
        omero.grid.DoubleColumn('Max', 'Max intensity', []),
        omero.grid.DoubleColumn('Sum', 'Integrated density', []),
        omero.grid.LongColumn('Area', 'Area in pixels', [])
    ]
    table.initialize(columns)

    # Add data
    data = [
        omero.grid.ImageColumn('Image', 'Source image', image_ids),
        omero.grid.LongColumn('ShapeID', 'ROI shape ID', shape_ids),
        omero.grid.DoubleColumn('Mean', 'Mean intensity', means),
        omero.grid.DoubleColumn('Min', 'Min intensity', mins),
        omero.grid.DoubleColumn('Max', 'Max intensity', maxs),
        omero.grid.DoubleColumn('Sum', 'Integrated density', sums),
        omero.grid.LongColumn('Area', 'Area in pixels', areas)
    ]
    table.addData(data)

    orig_file = table.getOriginalFile()
    table.close()

    # Link to dataset
    file_ann = omero.model.FileAnnotationI()
    file_ann.setFile(orig_file)
    file_ann = conn.getUpdateService().saveAndReturnObject(file_ann)

    link = omero.model.DatasetAnnotationLinkI()
    link.setParent(dataset._obj)
    link.setChild(file_ann)
    conn.getUpdateService().saveAndReturnObject(link)

    return file_ann

if __name__ == "__main__":
    run_script()
```

## 脚本部署

### 安装位置

脚本应放置在 OMERO 服务器脚本目录中：
```
OMERO_DIR/lib/scripts/
```

### 推荐结构

```
lib/scripts/
├── analysis/
│   ├── Cell_Counter.py
│   └── ROI_Analyzer.py
├── export/
│   ├── Export_Images.py
│   └── Export_ROIs.py
└── util/
    └── Helper_Functions.py
```

### 测试脚本

```bash
# Test script syntax
python Script_Name.py

# Upload to OMERO
omero script upload Script_Name.py

# List scripts
omero script list

# Run script from CLI
omero script launch Script_ID Dataset_ID=123
```

## 最佳实践

1. **错误处理**：始终使用try-finally来关闭会话
2. **进度更新**：打印长时间操作的状态消息
3. **参数验证**：处理前检查参数
4. **内存管理**：批量处理大型数据集
5. **文档**：包括清晰的描述和参数文档
6. **版本控制**：在脚本中包含版本号
7. **命名空间**：为输出使用适当的命名空间
8. **返回对象**：返回创建的对象以供客户端显示
9. **日志记录**：使用 print() 获取服务器日志
10. **测试**：使用各种输入组合进行测试

## 常见模式

### 进度报告

```python
total = len(images)
for idx, image in enumerate(images):
    print(f"Processing {idx + 1}/{total}: {image.getName()}")
    # Process image
```

### 错误集合

```python
errors = []
for image in images:
    try:
        process_image(image)
    except Exception as e:
        errors.append(f"{image.getName()}: {str(e)}")

if errors:
    message = "Completed with errors:\n" + "\n".join(errors)
else:
    message = "All images processed successfully"
```

### 资源清理

```python
try:
    # Script processing
    pass
finally:
    # Clean up temporary files
    if os.path.exists(temp_file):
        os.remove(temp_file)
    client.closeSession()
```