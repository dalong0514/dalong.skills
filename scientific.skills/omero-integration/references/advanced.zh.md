<!-- 此文件由机器翻译自 advanced.md -->

# 高级功能

本参考涵盖高级 OMERO 操作，包括权限、删除、文件集和管理任务。

## 删除对象

### 等待删除

```python
# Delete objects and wait for completion
project_ids = [1, 2, 3]
conn.deleteObjects("Project", project_ids, wait=True)
print("Deletion complete")

# Delete without waiting (asynchronous)
conn.deleteObjects("Dataset", [dataset_id], wait=False)
```

### 通过回调监听删除

<<<代码块_1>>>

### 删除不同的对象类型

<<<代码块_2>>>

### 使用级联删除

<<<代码块_3>>>

## 文件集

文件集表示原始导入文件的集合。它们是在 OMERO 5.0 中引入的。

### 检查图像是否有文件集

<<<代码块_4>>>

### 访问文件集信息

<<<代码块_5>>>

### 直接获取文件集

<<<代码块_6>>>

### 下载原始文件

```python
import os

fileset = image.getFileset()

if fileset:
    download_dir = "./original_files"
    os.makedirs(download_dir, exist_ok=True)

    for orig_file in fileset.listFiles():
        file_name = orig_file.getName()
        file_path = os.path.join(download_dir, file_name)

        print(f"Downloading: {file_name}")

        # Get file as RawFileStore
        raw_file_store = conn.createRawFileStore()
        raw_file_store.setFileId(orig_file.getId())

        # Download in chunks
        with open(file_path, 'wb') as f:
            offset = 0
            chunk_size = 1024 * 1024  # 1MB chunks
            size = orig_file.getSize()

            while offset < size:
                chunk = raw_file_store.read(offset, chunk_size)
                f.write(chunk)
                offset += len(chunk)

        raw_file_store.close()
        print(f"Saved to: {file_path}")
```

## 组权限

OMERO 使用基于组的权限来控制数据访问。

### 权限级别

- **PRIVATE** (`rw----`)：只有所有者可以读/写
- **READ-ONLY** (`rwr---`)：组成员可以读取，只有所有者可以写入
- **READ-ANNOTATE** (`rwra--`)：群组成员可以阅读和注释
- **READ-WRITE** (`rwrw--`)：组成员可以读写

### 检查当前组权限

```python
# Get current group
group = conn.getGroupFromContext()

# Get permissions
permissions = group.getDetails().getPermissions()
perm_string = str(permissions)

# Map to readable names
permission_names = {
    'rw----': 'PRIVATE',
    'rwr---': 'READ-ONLY',
    'rwra--': 'READ-ANNOTATE',
    'rwrw--': 'READ-WRITE'
}

perm_name = permission_names.get(perm_string, 'UNKNOWN')
print(f"Group: {group.getName()}")
print(f"Permissions: {perm_name} ({perm_string})")
```

### 列出用户组

```python
# Get all groups for current user
print("User's groups:")
for group in conn.getGroupsMemberOf():
    print(f"  {group.getName()} (ID: {group.getId()})")

    # Get group permissions
    perms = group.getDetails().getPermissions()
    print(f"    Permissions: {perms}")
```

### 获取群组成员

```python
# Get group
group = conn.getObject("ExperimenterGroup", group_id)

# List members
print(f"Members of {group.getName()}:")
for member in group.getMembers():
    print(f"  {member.getFullName()} ({member.getOmeName()})")
```

## 跨组查询

### 跨所有组查询

```python
# Set context to query all accessible groups
conn.SERVICE_OPTS.setOmeroGroup('-1')

# Now queries span all groups
image = conn.getObject("Image", image_id)
if image:
    group = image.getDetails().getGroup()
    print(f"Image found in group: {group.getName()}")

# List projects across all groups
for project in conn.getObjects("Project"):
    group = project.getDetails().getGroup()
    print(f"Project: {project.getName()} (Group: {group.getName()})")
```

### 切换到特定组

```python
# Get image's group
image = conn.getObject("Image", image_id)
group_id = image.getDetails().getGroup().getId()

# Switch to that group's context
conn.SERVICE_OPTS.setOmeroGroup(group_id)

# Subsequent operations use this group
projects = conn.listProjects()  # Only from this group
```

### 重置为默认组

```python
# Get default group
default_group_id = conn.getEventContext().groupId

# Switch back to default
conn.SERVICE_OPTS.setOmeroGroup(default_group_id)
```

## 行政运作

### 检查管理员状态

```python
# Check if current user is admin
if conn.isAdmin():
    print("User has admin privileges")

# Check if full admin
if conn.isFullAdmin():
    print("User is full administrator")
else:
    # Check specific privileges
    privileges = conn.getCurrentAdminPrivileges()
    print(f"Admin privileges: {privileges}")
```

### 列出管理员

```python
# Get all administrators
print("Administrators:")
for admin in conn.getAdministrators():
    print(f"  ID: {admin.getId()}")
    print(f"  Username: {admin.getOmeName()}")
    print(f"  Full Name: {admin.getFullName()}")
```

### 设置对象所有者（仅限管理员）

```python
import omero.model

# Create annotation with specific owner (requires admin)
tag_ann = omero.gateway.TagAnnotationWrapper(conn)
tag_ann.setValue("Admin-created tag")

# Set owner
user_id = 5
tag_ann._obj.details.owner = omero.model.ExperimenterI(user_id, False)
tag_ann.save()

print(f"Created annotation owned by user {user_id}")
```

### 替代用户连接（仅限管理员）

```python
# Connect as admin
admin_conn = BlitzGateway(admin_user, admin_pass, host=host, port=4064)
admin_conn.connect()

# Get target user
target_user_id = 10
user = admin_conn.getObject("Experimenter", target_user_id)
username = user.getOmeName()

# Create connection as that user
user_conn = admin_conn.suConn(username)

print(f"Connected as {username}")

# Perform operations as that user
for project in user_conn.listProjects():
    print(f"  {project.getName()}")

# Close connections
user_conn.close()
admin_conn.close()
```

### 列出所有用户

```python
# Get all users (admin operation)
print("All users:")
for user in conn.getObjects("Experimenter"):
    print(f"  ID: {user.getId()}")
    print(f"  Username: {user.getOmeName()}")
    print(f"  Full Name: {user.getFullName()}")
    print(f"  Email: {user.getEmail()}")
    print()
```

## 服务访问

OMERO 为特定操作提供各种服务。

### 更新服务

```python
# Get update service
updateService = conn.getUpdateService()

# Save and return object
roi = omero.model.RoiI()
roi.setImage(image._obj)
saved_roi = updateService.saveAndReturnObject(roi)

# Save multiple objects
objects = [obj1, obj2, obj3]
saved_objects = updateService.saveAndReturnArray(objects)
```

### 投资回报率服务

```python
# Get ROI service
roi_service = conn.getRoiService()

# Find ROIs for image
result = roi_service.findByImage(image_id, None)

# Get shape statistics
shape_ids = [shape.id.val for roi in result.rois
             for shape in roi.copyShapes()]
stats = roi_service.getShapeStatsRestricted(shape_ids, 0, 0, [0])
```

### 元数据服务

```python
# Get metadata service
metadataService = conn.getMetadataService()

# Load annotations by type and namespace
ns_to_include = ["mylab.analysis"]
ns_to_exclude = []

annotations = metadataService.loadSpecifiedAnnotations(
    'omero.model.FileAnnotation',
    ns_to_include,
    ns_to_exclude,
    None
)

for ann in annotations:
    print(f"Annotation: {ann.getId().getValue()}")
```

### 查询服务

```python
# Get query service
queryService = conn.getQueryService()

# Build query (more complex queries)
params = omero.sys.ParametersI()
params.addLong("image_id", image_id)

query = "select i from Image i where i.id = :image_id"
image = queryService.findByQuery(query, params)
```

### 缩略图服务

```python
# Get thumbnail service
thumbnailService = conn.createThumbnailStore()

# Set current image
thumbnailService.setPixelsId(image.getPrimaryPixels().getId())

# Get thumbnail
thumbnail = thumbnailService.getThumbnail(96, 96)

# Close service
thumbnailService.close()
```

### 原始文件存储

```python
# Get raw file store
rawFileStore = conn.createRawFileStore()

# Set file ID
rawFileStore.setFileId(orig_file_id)

# Read file
data = rawFileStore.read(0, rawFileStore.size())

# Close
rawFileStore.close()
```

## 对象所有权和详细信息

### 获取对象详细信息

```python
image = conn.getObject("Image", image_id)

# Get details
details = image.getDetails()

# Owner information
owner = details.getOwner()
print(f"Owner ID: {owner.getId()}")
print(f"Username: {owner.getOmeName()}")
print(f"Full Name: {owner.getFullName()}")

# Group information
group = details.getGroup()
print(f"Group: {group.getName()} (ID: {group.getId()})")

# Creation information
creation_event = details.getCreationEvent()
print(f"Created: {creation_event.getTime()}")

# Update information
update_event = details.getUpdateEvent()
print(f"Updated: {update_event.getTime()}")
```

### 获取权限

```python
# Get object permissions
details = image.getDetails()
permissions = details.getPermissions()

# Check specific permissions
can_edit = permissions.canEdit()
can_annotate = permissions.canAnnotate()
can_link = permissions.canLink()
can_delete = permissions.canDelete()

print(f"Can edit: {can_edit}")
print(f"Can annotate: {can_annotate}")
print(f"Can link: {can_link}")
print(f"Can delete: {can_delete}")
```

## 事件上下文

### 获取当前事件上下文

```python
# Get event context (current session info)
ctx = conn.getEventContext()

print(f"User ID: {ctx.userId}")
print(f"Username: {ctx.userName}")
print(f"Group ID: {ctx.groupId}")
print(f"Group Name: {ctx.groupName}")
print(f"Session ID: {ctx.sessionId}")
print(f"Is Admin: {ctx.isAdmin}")
```

## 完整的管理示例

```python
from omero.gateway import BlitzGateway

# Connect as admin
ADMIN_USER = 'root'
ADMIN_PASS = 'password'
HOST = 'omero.example.com'
PORT = 4064

with BlitzGateway(ADMIN_USER, ADMIN_PASS, host=HOST, port=PORT) as admin_conn:
    print("=== Administrator Operations ===\n")

    # List all users
    print("All Users:")
    for user in admin_conn.getObjects("Experimenter"):
        print(f"  {user.getOmeName()}: {user.getFullName()}")

    # List all groups
    print("\nAll Groups:")
    for group in admin_conn.getObjects("ExperimenterGroup"):
        perms = group.getDetails().getPermissions()
        print(f"  {group.getName()}: {perms}")

        # List members
        for member in group.getMembers():
            print(f"    - {member.getOmeName()}")

    # Query across all groups
    print("\nAll Projects (all groups):")
    admin_conn.SERVICE_OPTS.setOmeroGroup('-1')

    for project in admin_conn.getObjects("Project"):
        owner = project.getDetails().getOwner()
        group = project.getDetails().getGroup()
        print(f"  {project.getName()}")
        print(f"    Owner: {owner.getOmeName()}")
        print(f"    Group: {group.getName()}")

    # Connect as another user
    target_user_id = 5
    user = admin_conn.getObject("Experimenter", target_user_id)

    if user:
        print(f"\n=== Operating as {user.getOmeName()} ===\n")

        user_conn = admin_conn.suConn(user.getOmeName())

        # List that user's projects
        for project in user_conn.listProjects():
            print(f"  {project.getName()}")

        user_conn.close()
```

## 最佳实践

1. **权限**：操作前务必检查权限
2. **组上下文**：为查询设置适当的组上下文
3. **管理员操作**：谨慎谨慎地使用管理员权限
4. **删除确认**：删除对象之前始终确认
5. **回调监控**：通过回调监控长删除操作
6. **文件集感知**：处理图像时检查文件集
7. **服务清理**：完成后关闭服务（thumbnailStore、rawFileStore）
8. **跨组查询**：使用`-1`组ID进行跨组访问
9. **错误处理**：始终处理权限和访问错误
10. **文档**：清楚地记录管理操作

## 故障排除

### 权限被拒绝

```python
try:
    conn.deleteObjects("Project", [project_id], wait=True)
except Exception as e:
    if "SecurityViolation" in str(e):
        print("Permission denied: You don't own this object")
    else:
        raise
```

### 未找到对象

```python
# Check if object exists before accessing
obj = conn.getObject("Image", image_id)
if obj is None:
    print(f"Image {image_id} not found or not accessible")
else:
    # Process object
    pass
```

### 群体背景问题

```python
# If object not found, try cross-group query
conn.SERVICE_OPTS.setOmeroGroup('-1')
obj = conn.getObject("Image", image_id)

if obj:
    # Switch to object's group for further operations
    group_id = obj.getDetails().getGroup().getId()
    conn.SERVICE_OPTS.setOmeroGroup(group_id)
```