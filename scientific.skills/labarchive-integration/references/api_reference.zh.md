<!-- 此文件由机器翻译自 api_reference.md -->

# LabArchives API 参考

## API 结构

所有 LabArchives API 调用都遵循以下 URL 模式：

```
https://<base_url>/api/<api_class>/<api_method>?<authentication_parameters>&<method_parameters>
```

## 区域 API 端点

|地区 |基本网址 |
|--------|----------|
|美国/国际 | `https://api.labarchives.com/api` |
|澳大利亚 | `https://auapi.labarchives.com/api` |
|英国 | `https://ukapi.labarchives.com/api` |

## 身份验证

所有 API 调用都需要身份验证参数：

- `access_key_id`：由 LabArchives 管理员提供
- `access_password`：由 LabArchives 管理员提供
- 某些操作可能需要额外的用户特定凭据

## API 类和方法

### 用户 API 类

#### `users/user_access_info`

检索用户 ID 和笔记本访问信息。

**参数：**
- `login_or_email`（必填）：用户的电子邮件地址或登录用户名
- `password`（必需）：用户的外部应用程序密码（不是常规登录密码）

**返回：** XML 或 JSON 响应包含：
- 用户ID（uid）
- 具有 ID 的可访问笔记本列表 (nbid)
- 帐户状态和权限

**示例：**
<<<代码块_1>>>

#### `users/user_info_via_id`

通过用户ID检索详细的用户信息。

**参数：**
- `uid`（必需）：从 user_access_info 获取的用户 ID

**返回：** 用户个人资料信息包括：
- 姓名和电子邮件
- 帐户创建日期
- 机构隶属关系
- 角色和权限
- 存储配额和使用情况

**示例：**
<<<代码块_2>>>

### 笔记本 API 类

#### `notebooks/notebook_backup`

下载完整的笔记本数据，包括条目、附件和元数据。

**参数：**
- `uid`（必需）：用户 ID
- `nbid`（必需）：笔记本 ID
- `json`（可选，默认值：false）：以 JSON 格式而不是 XML 格式返回数据
- `no_attachments`（可选，默认值：false）：从备份中排除附件

**退货：**
- 当 `no_attachments=false` 时：包含所有笔记本数据的 7z 压缩存档
- 当`no_attachments=true`时：带有条目内容的XML或JSON结构化数据

**文件格式：**
返回的存档包括：
- HTML格式的输入文本内容
- 原始格式的文件附件
- 元数据 XML 文件，包含时间戳、作者和版本历史记录
- 评论主题和注释

**示例：**
<<<代码块_3>>>

<<<代码块_4>>>

#### `notebooks/list_notebooks`

检索用户可访问的所有笔记本（方法名称可能因 API 版本而异）。

**参数：**
- `uid`（必需）：用户 ID

**返回：** 笔记本列表，其中包含：
- 笔记本 ID (nbid)
- 笔记本名称
- 创建和修改日期
- 访问级别（所有者、编辑者、查看者）
- 会员数量

### 条目 API 类

#### `entries/create_entry`

在笔记本中创建一个新条目。

**参数：**
- `uid`（必需）：用户 ID
- `nbid`（必需）：笔记本 ID
- `title`（必填）：条目标题
- `content`（可选）：HTML 格式的条目内容
- `date`（可选）：输入日期（默认为当前日期）

**返回：** 条目 ID 和创建确认

**示例：**
<<<代码块_5>>>

#### `entries/create_comment`

向现有条目添加评论。

**参数：**
- `uid`（必需）：用户 ID
- `nbid`（必需）：笔记本 ID
- `entry_id`（必需）：目标条目 ID
- `comment`（必需）：注释文本（支持 HTML）

**返回：** 评论 ID 和时间戳

#### `entries/create_part`

将组件/部分添加到条目（例如文本部分、表格、图像）。

**参数：**
- `uid`（必需）：用户 ID
- `nbid`（必需）：笔记本 ID
- `entry_id`（必需）：目标条目 ID
- `part_type`（必需）：部分类型（文本、表格、图像等）
- `content`（必需）：采用适当格式的部分内容

**返回：** 零件 ID 和创建确认

#### `entries/upload_attachment`

将文件附件上传到条目。

**参数：**
- `uid`（必需）：用户 ID
- `nbid`（必需）：笔记本 ID
- `entry_id`（必需）：目标条目 ID
- `file`（必需）：文件数据（多部分/表单数据）
- `filename`（必需）：原始文件名

**返回：** 附件 ID 和上传确认

**使用请求库的示例：**
<<<代码块_6>>>

### 站点报告 API 类

仅限企业的机构报告和分析功能。

#### `site_reports/detailed_usage_report`

为机构生成全面的使用统计数据。
**参数：**
- `start_date`（必需）：报告开始日期 (YYYY-MM-DD)
- `end_date`（必需）：报告结束日期 (YYYY-MM-DD)
- `format`（可选）：输出格式（csv、json、xml）

**回报：** 使用指标包括：
- 用户登录频率
- 条目创建计数
- 存储利用率
- 协作统计
- 基于时间的活动模式

#### `site_reports/detailed_notebook_report`

生成机构内所有笔记本的详细报告。

**参数：**
- `include_settings`（可选，默认值：false）：包括笔记本设置
- `include_members`（可选，默认值：false）：包括成员列表

**退货：** 笔记本库存包含：
- 笔记本名称和 ID
- 业主信息
- 创建和最后修改日期
- 会员数量和访问级别
- 存储大小
- 设置（如果需要）

#### `site_reports/pdf_offline_generation_report`

出于合规性和审核目的跟踪 PDF 导出。

**参数：**
- `start_date`（必需）：报告开始日期
- `end_date`（必需）：报告结束日期

**返回：** 导出活动日志：
- 生成PDF的用户
- 笔记本和条目导出
- 导出时间戳
- IP地址

### 实用程序 API 类

#### `utilities/institutional_login_urls`

检索机构登录 URL 以进行 SSO 集成。

**参数：** 不需要（使用访问密钥身份验证）

**返回：** 机构登录端点列表

## 响应格式

### XML 响应示例

```xml
<?xml version="1.0" encoding="UTF-8"?>
<response>
    <uid>12345</uid>
    <email>researcher@university.edu</email>
    <notebooks>
        <notebook>
            <nbid>67890</nbid>
            <name>Lab Notebook 2025</name>
            <role>owner</role>
        </notebook>
    </notebooks>
</response>
```

### JSON 响应示例

```json
{
    "uid": "12345",
    "email": "researcher@university.edu",
    "notebooks": [
        {
            "nbid": "67890",
            "name": "Lab Notebook 2025",
            "role": "owner"
        }
    ]
}
```

## 错误代码

|代码|留言 |意义|解决方案 |
|------|---------|---------|---------|
| 401 | 401未经授权 |凭证无效 |验证 access_key_id 和 access_password |
| 403 | 403禁止 |权限不足|检查用户角色和笔记本访问权限 |
| 404 | 404未找到 |资源不存在 |验证 uid、nbid 或entry_id 是否正确 |
| 429 | 429太多请求 |超出速率限制 |实施指数退避 |
| 500 | 500内部服务器错误 |服务器端问题 |重试请求或联系支持人员 |

## 速率限制

LabArchives 实施速率限制以保证服务稳定性：

- **推荐：** 每个 API 密钥每分钟最多 60 个请求
- **突发限额：** 可以容忍最多 100 个请求的短突发
- **最佳实践：** 在批处理操作的请求之间实现 1-2 秒的延迟

## API 版本控制

LabArchives API 向后兼容。添加新方法不会破坏现有的实现。监视 LabArchives 公告以了解新功能。

## 支持和文档

对于 API 访问请求、技术问题或功能请求：
- 电子邮件：support@labarchives.com
- 包括您的机构名称和具体用例以获得更快的帮助