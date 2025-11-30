<!-- 此文件由机器翻译自 integrations.md -->

# LabArchives 第三方集成

## 概述

LabArchives 与众多科学软件平台集成，以简化研究工作流程。本文档涵盖了每个受支持平台的编程集成方法、自动化策略和最佳实践。

## 集成类别

### 1. 协议管理

#### Protocols.io 集成

将协议直接从 Protocols.io 导出到 LabArchives 笔记本。

**使用案例：**
- 标准化实验室笔记本的实验程序
- 维护协议的版本控制
- 将协议与实验结果链接起来

**设置：**
1. 在 LabArchives 设置中启用 Protocols.io 集成
2. 使用 Protocols.io 帐户进行身份验证
3. 浏览并选择要导出的协议

**程序化方法：**
```python
# Export Protocols.io protocol as HTML/PDF
# Then upload to LabArchives via API

def import_protocol_to_labarchives(client, uid, nbid, protocol_id):
    """Import Protocols.io protocol to LabArchives entry"""
    # 1. Fetch protocol from Protocols.io API
    protocol_data = fetch_protocol_from_protocolsio(protocol_id)

    # 2. Create new entry in LabArchives
    entry_params = {
        'uid': uid,
        'nbid': nbid,
        'title': f"Protocol: {protocol_data['title']}",
        'content': protocol_data['html_content']
    }
    response = client.make_call('entries', 'create_entry', params=entry_params)

    # 3. Add protocol metadata as comment
    entry_id = extract_entry_id(response)
    comment_params = {
        'uid': uid,
        'nbid': nbid,
        'entry_id': entry_id,
        'comment': f"Protocols.io ID: {protocol_id}<br>Version: {protocol_data['version']}"
    }
    client.make_call('entries', 'create_comment', params=comment_params)

    return entry_id
```

**更新日期：** 2025 年 9 月 22 日

### 2. 数据分析工具

#### GraphPad Prism 集成（版本 8+）

将分析、图表和数据直接从 Prism 导出到 LabArchives。

**使用案例：**
- 使用原始数据存档统计分析
- 出版物的文档图形生成
- 维护分析审计跟踪以确保合规性

**设置：**
1.安装GraphPad Prism 8或更高版本
2. 在 Prism 首选项中配置 LabArchives 连接
3. 使用“文件”菜单中的“导出到 LabArchives”选项

**程序化方法：**
<<<代码块_1>>>

**支持的文件类型：**
- .pzfx（棱镜项目文件）
- .png、.jpg、.pdf（导出的图表）
- .xlsx（导出的数据表）

**更新日期：** 2025 年 9 月 8 日

### 3.分子生物学与生物信息学

#### SnapGene 集成

直接集成分子生物学工作流程、质粒图谱和序列分析。

**使用案例：**
- 记录克隆策略
- 存档质粒图谱和实验记录
- 将序列与实验结果链接起来

**设置：**
1.安装SnapGene软件
2. 在 SnapGene 首选项中启用 LabArchives 导出
3.使用“发送到LabArchives”功能

**文件格式支持：**
- .dna（SnapGene 文件）
- .gb、.gbk（GenBank 格式）
- .fasta（序列文件）
- .png、.pdf（质粒图导出）

**程序化工作流程：**
<<<代码块_2>>>

#### 慷慨的整合

生物信息学分析从 Geneious 导出到 LabArchives。

**使用案例：**
- 存档序列比对和系统发育树
- 记录 NGS 分析流程
- 将生物信息学工作流程与湿实验室实验联系起来

**支持导出：**
- 序列比对
- 系统发育树
- 大会报告
- 变异调用结果

**文件格式：**
- .geneious（慷慨的文件）
- .fasta、.fastq（序列数据）
- .bam、.sam（对齐文件）
- .vcf（变体文件）

### 4. 计算笔记本

#### Jupyter 集成

将 Jupyter 笔记本嵌入为 LabArchives 条目，以进行可重复的计算研究。

**使用案例：**
- 记录数据分析工作流程
- 存档计算实验
- 链接代码、结果和叙述

**工作流程：**

<<<代码块_3>>>

**最佳实践：**
- 包含输出的导出（导出前运行所有单元）
- 包含environment.yml或requirements.txt作为附件
- 在注释中添加执行时间戳和系统信息

### 5. 临床研究

#### REDCap 集成

临床数据采集与 LabArchives 集成，以实现研究合规性和审计追踪。

**使用案例：**
- 将临床数据收集链接到研究笔记本
- 维护监管合规性的审计跟踪
- 记录临床试验方案和修订

**整合方式：**
- REDCap API 将数据导出到 LabArchives 条目
- 纵向研究的自动数据同步
- 符合 HIPAA 的数据处理

**工作流程示例：**
<<<代码块_4>>>

**合规特性：**
- 21 CFR 第 11 部分合规性
- 审计跟踪维护
- 数据完整性验证

### 6. 研究出版

#### Qeios 集成

用于预印本和同行评审的研究出版平台集成。

**使用案例：**
- 将研究成果导出到预印本服务器
- 文档发布工作流程
- 将已发表的文章链接到实验室笔记本

**工作流程：**
- 从 LabArchives 导出格式化条目
- 提交至Qeios平台
- 维护笔记本和出版物之间的双向链接

#### SciSpace 集成

文献管理和引文整合。

**使用案例：**
- 链接参考实验程序
- 在笔记本中维护文献综述
- 生成报告参考书目

**特点：**
- 引文从 SciSpace 导入到 LabArchives
- PDF注释同步
- 参考管理
## 集成的 OAuth 身份验证

LabArchives 现在使用 OAuth 2.0 进行新的第三方集成。

**应用程序开发人员的 OAuth 流程：**

<<<代码块_5>>>

**OAuth 优点：**
- 比 API 密钥更安全
- 细粒度的权限控制
- 长期运行集成的令牌刷新
- 可撤销的访问权限

## 定制集成开发

### 一般工作流程

对于未正式支持的工具，开发自定义集成：

1. **从源应用程序导出数据**（API或文件导出）
2. **将格式**转换为 HTML 或支持的文件类型
3. **使用 LabArchives API 进行身份验证**
4. **创建条目**或上传附件
5. **通过注释添加元数据**以实现可追溯性

### 示例：自定义集成模板

<<<代码块_6>>>

## 集成最佳实践

1. **版本控制：** 跟踪哪个软件版本生成了数据
2. **元数据保存：** 包括时间戳、用户信息和处理参数
3. **文件格式标准：** 尽可能使用开放格式（CSV、JSON、HTML）
4. **批量操作：** 对批量上传实施限速
5. **错误处理：** 使用指数退避实现重试逻辑
6. **审核跟踪：** 记录所有 API 操作以确保合规性
7. **测试：** 在生产使用之前验证测试笔记本中的集成

## 集成故障排除

### 常见问题

**集成未出现在 LabArchives 中：**
- 验证管理员是否启用了集成
- 如果使用 OAuth，请检查 OAuth 权限
- 确保兼容的软件版本

**文件上传失败：**
- 验证文件大小限制（通常每个文件 2GB）
- 检查文件格式兼容性
- 确保足够的存储配额

**身份验证错误：**
- 验证 API 凭证是否最新
- 检查集成特定令牌是否已过期
- 确认用户拥有必要的权限

### 集成支持

对于特定于集成的问题：
- 检查软件供应商文档（例如 GraphPad、Protocols.io）
- 联系 LabArchives 支持：support@labarchives.com
- 查看 LabArchives 知识库：help.labarchives.com

## 未来的整合机会

定制开发的潜在集成：
- 电子数据采集（EDC）系统
- 实验室信息管理系统（LIMS）
- 仪器数据系统（色谱、光谱）
- 云存储平台（Box、Dropbox、Google Drive）
- 项目管理工具（Asana、Monday.com）
- 赠款管理系统

对于定制集成开发，请联系 LabArchives 获取 API 合作机会。