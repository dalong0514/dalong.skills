<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：实验室档案集成
描述：“电子实验室笔记本 API 集成。访问笔记本、管理条目/附件、备份笔记本、与 Protocols.io/Jupyter/REDCap 集成，实现编程式 ELN 工作流程。”
---

# 实验室档案集成

## 概述

LabArchives 是一个用于研究文档和数据管理的电子实验室笔记本平台。访问笔记本、管理条目和附件、生成报告以及通过 REST API 以编程方式与第三方工具集成。

## 何时使用此技能

该技能应该在以下情况下使用：
- 使用 LabArchives REST API 实现笔记本自动化
- 以编程方式备份笔记本
- 创建或管理笔记本条目和附件
- 生成站点报告和分析
- 将 LabArchives 与第三方工具集成（Protocols.io、Jupyter、REDCap）
- 自动将数据上传到电子实验室笔记本
- 以编程方式管理用户访问和权限

## 核心能力

### 1. 身份验证和配置

为 LabArchives API 集成设置 API 访问凭据和区域端点。

**先决条件：**
- 启用 API 访问的 Enterprise LabArchives 许可证
- LabArchives 管理员提供的 API 访问密钥 ID 和密码
- 用户身份验证凭据（电子邮件和外部应用程序密码）

**配置设置：**

使用 `scripts/setup_config.py` 脚本创建配置文件：

```bash
python3 scripts/setup_config.py
```

这将创建一个具有以下结构的 `config.yaml` 文件：

<<<代码块_1>>>

**区域 API 端点：**
- 美国/国际：`https://api.labarchives.com/api`
- 澳大利亚：`https://auapi.labarchives.com/api`
- 英国：`https://ukapi.labarchives.com/api`

有关详细的身份验证说明和故障排除，请参阅`references/authentication_guide.md`。

### 2.用户信息检索

获取后续API操作所需的用户ID（UID）和访问信息。

**工作流程：**

1. 使用登录凭据调用 `users/user_access_info` API 方法
2. 解析 XML/JSON 响应以提取用户 ID (UID)
3.使用UID通过`users/user_info_via_id`检索详细的用户信息

**使用 Python 包装器的示例：**

<<<代码块_2>>>

### 3.笔记本操作

管理笔记本访问、备份和元数据检索。

**关键操作：**

- **列出笔记本：** 检索用户可访问的所有笔记本
- **备份笔记本：** 下载完整的笔记本数据以及可选的附件包含
- **获取笔记本 ID：** 检索机构定义的笔记本标识符以与赠款/项目管理系统集成
- **获取笔记本成员：** 列出有权访问特定笔记本的所有用户
- **获取笔记本设置：** 检索笔记本的配置和权限

**笔记本备份示例：**

使用 `scripts/notebook_operations.py` 脚本：

<<<代码块_3>>>

**API端点格式：**
<<<代码块_4>>>

有关全面的 API 方法文档，请参阅 `references/api_reference.md`。

### 4. 条目和附件管理

创建、修改和管理笔记本条目和文件附件。

**录入操作：**
- 在笔记本中创建新条目
- 为现有条目添加评论
- 创建条目零件/组件
- 上传文件附件到条目

**附件工作流程：**

使用 `scripts/entry_operations.py` 脚本：

<<<代码块_5>>>

**支持的文件类型：**
- 文档（PDF、DOCX、TXT）
- 图片（PNG、JPG、TIFF）
- 数据文件（CSV、XLSX、HDF5）
- 科学格式（CIF、MOL、PDB）
- 档案（ZIP、7Z）

### 5. 站点报告和分析

生成有关笔记本使用情况、活动和合规性的机构报告（企业功能）。

**可用报告：**
- 详细的使用情况报告：用户活动指标和参与度统计数据
- 详细的笔记本报告：笔记本元数据、成员列表和设置
- PDF/离线笔记本生成报告：导出跟踪以确保合规性
- 笔记本成员报告：访问控制和协作分析
- 笔记本设置报告：配置和权限审核

**报告生成：**

<<<代码块_6>>>

### 6. 第三方集成

LabArchives 与众多科学软件平台集成。此技能提供了以编程方式利用这些集成的指导。

**支持的集成：**
- **Protocols.io：** 将协议直接导出到 LabArchives 笔记本
- **GraphPad Prism：** 导出分析和数据（版本 8+）
- **SnapGene：** 直接分子生物学工作流程集成
- **Geneious：** 生物信息学分析导出
- **Jupyter：** 嵌入 Jupyter 笔记本作为条目
- **REDCap：** 临床数据采集集成
- **Qeios：** 研究发布平台
- **SciSpace：** 文献管理

**OAuth 身份验证：**
LabArchives 现在使用 OAuth 进行所有新集成。旧集成可能使用 API 密钥身份验证。

有关详细的集成设置说明和用例，请参阅`references/integrations.md`。

## 常见工作流程

### 完整的笔记本备份工作流程

1. 认证并获取用户ID
2. 列出所有可访问的笔记本
3. 遍历笔记本并备份每一个
4. 使用时间戳元数据存储备份

```bash
# Complete backup script
python3 scripts/notebook_operations.py backup-all --email user@example.edu --password AUTH_TOKEN
```

### 自动数据上传工作流程

1. 使用 LabArchives API 进行身份验证
2. 识别目标笔记本和条目
3.上传实验数据文件
4.为条目添加元数据注释
5. 生成活动报告

### 集成工作流程示例（Jupyter → LabArchives）

1. 将 Jupyter Notebook 导出为 HTML 或 PDF
2.使用entry_operations.py上传到LabArchives
3. 添加带有执行时间戳和环境信息的注释
4.标签输入，方便检索

## Python 包安装

安装 `labarchives-py` 包装器以简化 API 访问：

```bash
git clone https://github.com/mcmero/labarchives-py
cd labarchives-py
uv pip install .
```

或者，通过 Python 的 `requests` 库使用直接 HTTP 请求进行自定义实现。

## 最佳实践

1. **速率限制：** 在 API 调用之间实施适当的延迟以避免限流
2. **错误处理：** 始终将 API 调用包装在带有适当日志记录的 try- except 块中
3. **身份验证安全性：** 将凭据存储在环境变量或安全配置文件中（切勿在代码中）
4. **备份验证：**笔记本备份后，验证文件的完整性和完整性
5. **增量操作：** 对于大型笔记本，使用分页和批处理
6. **区域端点：** 使用正确的区域 API 端点以获得最佳性能

## 故障排除

**常见问题：**

- **401 Unauthorized:** 验证访问密钥 ID 和密码是否正确；检查您的帐户是否启用了 API 访问
- **404 Not Found:** 确认笔记本 ID (nbid) 存在并且用户具有访问权限
- **403 Forbidden:** 检查请求操作的用户权限
- **空响应：** 确保正确提供所需参数（uid、nbid）
- **附件上传失败：**验证文件大小限制和格式兼容性

如需更多支持，请联系 LabArchives：support@labarchives.com。

## 资源

该技能包括支持 LabArchives API 集成的捆绑资源：

### 脚本/

- `setup_config.py`：API 凭证的交互式配置文件生成器
- `notebook_operations.py`：用于列出、备份和管理笔记本的实用程序
- `entry_operations.py`：用于创建条目和上传附件的工具

###参考资料/

- `api_reference.md`：包含参数和示例的综合 API 端点文档
- `authentication_guide.md`：详细的身份验证设置和配置说明
- `integrations.md`：第三方集成设置指南和用例