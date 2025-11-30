<!-- 此文件由机器翻译自 microscopy_imaging_formats.md -->

# 显微镜和成像文件格式参考

本参考涵盖显微镜、医学成像、遥感和科学图像分析中使用的文件格式。

## 显微镜专用格式

### .tif / .tiff - 标记图像文件格式
**描述：** 支持多页面和元数据的灵活图像格式
**典型数据：** 显微图像、z 堆栈、时间序列、多通道
**用例：** 荧光显微镜、共焦成像、生物成像
**Python 库：**
- `tifffile`: `tifffile.imread('file.tif')` - 显微镜 TIFF 支持
- `PIL/Pillow`: `Image.open('file.tif')` - 基本 TIFF
- `scikit-image`：`io.imread('file.tif')`
- `AICSImageIO`：多格式显微镜阅读器
**EDA 方法：**
- 图像尺寸和位深度
- 多页面/z 堆栈分析
- 元数据提取（OME-TIFF）
- 频道分析和强度分布
- 时间动态（延时）
- 像素大小和空间校准
- 每个通道的直方图分析
- 动态范围利用

### .nd2 - 尼康 NIS-Elements
**描述：** 尼康专有显微镜格式
**典型数据：** 多维显微镜 (XYZCT)
**使用案例：** 尼康显微镜数据、共焦、宽视野
**Python 库：**
- `nd2reader`: `ND2Reader('file.nd2')`
- `pims`：`pims.ND2_Reader('file.nd2')`
- `AICSImageIO`：通用阅读器
**EDA 方法：**
- 实验元数据提取
- 通道配置
- 延时帧分析
- Z 轴堆叠深度和间距
- XY 载物台位置
- 激光设置和功率
- 像素合并信息
- 获取时间戳

### .lif - 徕卡图像格式
**描述：** 徕卡显微镜专有格式
**典型数据：** 多实验、多维图像
**用例：** 徕卡共焦和宽场数据
**Python 库：**
- `readlif`：`readlif.LifFile('file.lif')`
- `AICSImageIO`：LIF 支持
- `python-bioformats`：通过生物格式
**EDA 方法：**
- 多次实验检测
- 图像系列枚举
- 每个实验的元数据
- 频道和时间点结构
- 物理尺寸提取
- 物镜和探测器信息
- 扫描设置分析

### .czi - 卡尔蔡司图像
**描述：**蔡司显微镜格式
**典型数据：** 具有丰富元数据的多维显微镜
**使用案例：**蔡司共焦、光片、宽场
**Python 库：**
- `czifile`：`czifile.CziFile('file.czi')`
- `AICSImageIO`：CZI 支持
- `pylibCZIrw`：官方蔡司库
**EDA 方法：**
- 场景和位置分析
- 马赛克瓷砖结构
- 通道波长信息
- 采集模式检测
- 缩放和校准
- 仪器配置
- 投资回报率定义

### .oib / .oif - 奥林巴斯图像格式
**描述：** 奥林巴斯显微镜格式
**典型数据：** 共焦和多光子成像
**使用案例：** Olympus FluoView 数据
**Python 库：**
- `AICSImageIO`：OIB/OIF 支持
- `python-bioformats`：通过生物格式
**EDA 方法：**
- 目录结构验证（OIF）
- 元数据文件解析
- 通道配置
- 扫描参数
- 目标和过滤器信息
- PMT设置

### .vsi - 奥林巴斯 VSI
**描述：** 奥林巴斯幻灯片扫描仪格式
**典型数据：** 整个载玻片成像，大马赛克
**用例：** 虚拟显微镜、病理学
**Python 库：**
- `openslide-python`：`openslide.OpenSlide('file.vsi')`
- `AICSImageIO`：VSI 支持
**EDA 方法：**
- 金字塔层次分析
- 瓷砖结构和重叠
- 宏和标签图像
- 放大倍数
- 整张幻灯片统计
- 区域检测

### .ims - Imaris 格式
**描述：** 基于位平面 Imaris HDF5 的格式
**典型数据：** 大型 3D/4D 显微镜数据集
**用例：** 3D 渲染、延时分析
**Python 库：**
- `h5py`：直接 HDF5 访问
- `imaris_ims_file_reader`：专业阅读器
**EDA 方法：**
- 分辨率级别分析
- 时间点结构
- 渠道组织
- 数据集层次结构
- 缩略图生成
- 内存映射访问策略
- 分块优化

### .lsm - 蔡司 LSM
**描述：** 传统蔡司共焦格式
**典型数据：** 共焦激光扫描显微镜
**用例：** 较旧的蔡司共焦数据
**Python 库：**
- `tifffile`：LSM 支持（基于 TIFF）
- `python-bioformats`：LSM 读取
**EDA 方法：**
- 类似于具有 LSM 特定元数据的 TIFF
- 扫描速度和分辨率
- 激光线和功率
- 探测器增益和偏移
- LUT信息
### .stk - MetaMorph 堆栈
**描述：** MetaMorph 图像堆栈格式
**典型数据：** 延时或 z 堆栈序列
**用例：** MetaMorph 软件输出
**Python 库：**
- `tifffile`：STK 基于 TIFF
- `python-bioformats`：STK 支持
**EDA 方法：**
- 堆栈维数
- 平面元数据
- 计时信息
- 舞台位置
- UIC标签解析

### .dv - DeltaVision
**描述：** 应用 Precision DeltaVision 格式
**典型数据：** 反卷积显微镜
**使用案例：** DeltaVision 显微镜数据
**Python 库：**
- `mrc`：可以读取 DV（MRC 相关）
- `AICSImageIO`：DV 支持
**EDA 方法：**
- 波形信息（通道）
- 扩展标头分析
- 镜头和放大倍率
- 反卷积状态
- 每个部分的时间戳

### .mrc - 医学研究委员会
**描述：** 电子显微镜格式
**典型数据：** EM 图像、冷冻电镜、断层扫描
**用例：**结构生物学、电子显微镜
**Python 库：**
- `mrcfile`：`mrcfile.open('file.mrc')`
- `EMAN2`：EM 专用工具
**EDA 方法：**
- 体积尺寸
- 体素大小和单位
- 原点和地图统计
- 对称信息
- 扩展标头分析
- 密度统计
- 标头一致性验证

### .dm3 / .dm4 - Gatan 数码显微照片
**描述：** Gatan TEM/STEM 格式
**典型数据：** 透射电子显微镜
**使用案例：** TEM 成像和分析
**Python 库：**
- `hyperspy`：`hs.load('file.dm3')`
- `ncempy`：`ncempy.io.dm.dmReader('file.dm3')`
**EDA 方法：**
- 显微镜参数
- 能量色散光谱数据
- 衍射图案
- 校准信息
- 标签结构分析
- 图像系列处理

### .eer - 电子事件表示
**描述：** 直接电子探测器格式
**典型数据：** 来自探测器的电子计数数据
**用例：** 冷冻电镜数据收集
**Python 库：**
- `mrcfile`：一些 EER 支持
- 供应商特定工具（Gatan、TFS）
**EDA 方法：**
- 事件计数统计
- 帧速率和剂量
- 探测器配置
- 运动矫正评估
- 获得参考验证

### .ser - TIA 系列
**描述：** FEI/TFS TIA 格式
**典型数据：** EM 图像系列
**使用案例：** FEI/Thermo Fisher EM 数据
**Python 库：**
- `hyperspy`：SER 支持
- `ncempy`：TIA 读取器
**EDA 方法：**
- 系列结构
- 校准数据
- 获取元数据
- 时间戳
- 多维度数据组织

## 医学和生物成像

### .dcm - DICOM
**描述：** 医学中的数字成像和通信
**典型数据：** 带有患者/研究元数据的医学图像
**使用案例：** 临床成像、放射学、CT、MRI、PET
**Python 库：**
- `pydicom`：`pydicom.dcmread('file.dcm')`
- `SimpleITK`：`sitk.ReadImage('file.dcm')`
- `nibabel`：有限的 DICOM 支持
**EDA 方法：**
- 患者元数据提取（匿名检查）
- 特定模态分析
- 系列和研究组织
- 切片厚度和间距
- 窗口/级别设置
- 亨斯菲尔德单位（CT）
- 图像方向和位置
- 多帧分析

### .nii / .nii.gz - NIfTI
**描述：** 神经影像信息学技术倡议
**典型数据：** 脑成像、功能磁共振成像、结构磁共振成像
**用例：** 神经影像研究、大脑分析
**Python 库：**
- `nibabel`：`nibabel.load('file.nii')`
- `nilearn`：使用 ML 进行神经影像学
- `SimpleITK`：NIfTI 支持
**EDA 方法：**
- 体积尺寸和体素尺寸
- 仿射变换矩阵
- 时间序列分析（fMRI）
- 强度分布
- 大脑提取质量
- 注册评估
- 方向验证
- 标头信息一致性

### .mnc - MINC 格式
**描述：** 医学图像 NetCDF
**典型数据：** 医学成像（NIfTI 的前身）
**用例：** 遗留神经影像数据
**Python 库：**
- `pyminc`：MINC 特定工具
- `nibabel`：MINC 支持
**EDA 方法：**
- 类似于 NIfTI
- NetCDF结构探索
- 尺寸排序
- 元数据提取

### .nrrd - 近乎原始的栅格数据
**描述：** 具有独立标题的医学成像格式
**典型数据：** 医学图像、研究成像
**用例：** 3D 切片器、基于 ITK 的应用程序
**Python 库：**
- `pynrrd`：`nrrd.read('file.nrrd')`
- `SimpleITK`：NRRD 支持
**EDA 方法：**
- 标头字段分析
- 编码格式
- 尺寸和间距
- 方向矩阵
- 压缩评估
- 字节序处理

### .mha / .mhd - 元图像
**描述：**元图像格式（ITK）
**典型数据：** 医学/科学 3D 图像
**用例：** ITK/SimpleITK 应用程序
**Python 库：**
- `SimpleITK`：本机 MHA/MHD 支持
- `itk`：直接 ITK 集成
**EDA 方法：**
- 标头数据文件配对 (MHD)
- 变换矩阵
- 元素间距
- 压缩格式
- 数据类型和维度

### .hdr / .img - 分析格式
**描述：** 传统医学成像格式
**典型数据：** 脑成像（NIfTI 之前）
**用例：**旧的神经影像数据集
**Python 库：**
- `nibabel`：分析支持
- 建议转换为 NIfTI
**EDA 方法：**
- 标题-图像配对验证
- 字节顺序问题
- 转换为现代格式
- 元数据限制

## 科学图像格式

### .png - 便携式网络图形
**描述：** 无损压缩图像格式
**典型数据：** 2D 图像、屏幕截图、处理后的数据
**用例：**出版数字、无损存储
**Python 库：**
- `PIL/Pillow`：`Image.open('file.png')`
- `scikit-image`：`io.imread('file.png')`
- `imageio`：`imageio.imread('file.png')`
**EDA 方法：**
- 位深度分析（8位、16位）
- 颜色模式（灰度、RGB、调色板）
- 元数据（PNG 块）
- 透明度处理
- 压缩效率
- 直方图分析

### .jpg / .jpeg - 联合摄影专家组
**描述：** 有损压缩图像格式
**典型数据：** 自然图像、照片
**用例：** 可视化、网页图形（非原始数据）
**Python 库：**
- `PIL/Pillow`：标准 JPEG 支持
- `scikit-image`：JPEG 读取
**EDA 方法：**
- 压缩伪影检测
- 品质因数估计
- 色彩空间（RGB、灰度）
- EXIF 元数据
- 量化表分析
- 注意：不适合定量分析

### .bmp - 位图图像
**描述：**未压缩的光栅图像
**典型数据：** 简单图像、屏幕截图
**使用案例：**兼容性、简单存储
**Python 库：**
- `PIL/Pillow`：BMP 支持
- `scikit-image`：BMP 读取
**EDA 方法：**
- 颜色深度
- 调色板分析（如果已索引）
- 文件大小效率
- 像素格式验证

### .gif - 图形交换格式
**描述：** 支持动画的图像格式
**典型数据：** 动画图像、简单图形
**用例：** 动画、延时可视化
**Python 库：**
- `PIL/Pillow`：GIF 支持
- `imageio`：更好的 GIF 动画支持
**EDA 方法：**
- 帧计数和计时
- 调色板限制（256 色）
- 循环计数
- 废弃处理方法
- 透明度处理

### .svg - 可缩放矢量图形
**描述：** 基于 XML 的矢量图形
**典型数据：** 矢量绘图、绘图、图表
**用例：** 出版质量的图表、图表
**Python 库：**
- `svgpathtools`：路径操作
- `cairosvg`：光栅化
- `lxml`：XML解析
**EDA 方法：**
- 元素结构分析
- 款式信息
- 视图框和尺寸
- 路径复杂度
- 文本元素提取
- 分层组织

### .eps - 封装的 PostScript
**说明：**矢量图形格式
**典型数据：** 出版物数据
**用例：** 旧出版物图形
**Python 库：**
- `PIL/Pillow`：基本 EPS 光栅化
- `ghostscript` 通过子进程
**EDA 方法：**
- 边界框信息
- 预览图像验证
- 字体嵌入
- 转换为现代格式

### .pdf（图片）
**描述：** 带图像的便携式文档格式
**典型数据：** 出版物数据、多页文档
**用例：** 出版、数据呈现
**Python 库：**
- `PyMuPDF/fitz`：`fitz.open('file.pdf')`
- `pdf2image`：光栅化
- `pdfplumber`：文本和布局提取
**EDA 方法：**
- 页数
- 图像提取
- 分辨率和DPI
- 嵌入字体和元数据
- 压缩方法
- 图像与矢量内容

### .fig - MATLAB 图
**说明：** MATLAB 图形文件
**典型数据：** MATLAB 绘图和图形
**用例：** MATLAB 数据可视化
**Python 库：**
- 自定义解析器（MAT 文件结构）
- 转换为其他格式
**EDA 方法：**
- 图形结构
- 从图中提取数据
- 轴和标签信息
- 地块类型识别

### .hdf5（特定于成像）
**描述：** HDF5 适用于大型成像数据集
**典型数据：** 高内涵筛选，大型显微镜
**用例：** BigDataViewer，大规模成像
**Python 库：**
- `h5py`：通用 HDF5 访问
- 成像专用阅读器 (BigDataViewer)
**EDA 方法：**
- 数据集层次结构
- 块和压缩策略
- 多分辨率金字塔
- 元数据组织
- 内存映射访问
- 并行I/O性能

### .zarr - 分块数组存储
**描述：** 云优化的阵列存储
**典型数据：** 大型成像数据集，OME-ZARR
**用例：** 云显微镜、大规模分析
**Python 库：**
- `zarr`：`zarr.open('file.zarr')`
- `ome-zarr-py`：OME-ZARR 支持
**EDA 方法：**
- 块大小优化
- 压缩编解码器分析
- 多尺度表示
- 数组维度和数据类型
- 元数据结构（OME）
- 云访问模式

### .raw - 原始图像数据
**描述：**未格式化的二进制像素数据
**典型数据：** 原始检测器输出
**用例：** 定制成像系统
**Python 库：**
- `numpy`: `np.fromfile()` 与数据类型
- `imageio`：原始格式插件
**EDA 方法：**
- 尺寸确定（需要外部信息）
- 字节顺序和数据类型
- 标头存在检测
- 像素值范围
- 噪声特性

### .bin - 二进制图像数据
**描述：** 通用二进制图像格式
**典型数据：** 原始或自定义格式的图像
**用例：** 仪器特定输出
**Python 库：**
- `numpy`：自定义二进制读取
- `struct`：用于结构化二进制数据
**EDA 方法：**
- 需要格式规范
- 标头解析（如果存在）
- 数据类型推断
- 维度提取
- 使用已知参数进行验证

## 图像分析格式

### .roi - ImageJ 投资回报率
**描述：** ImageJ 感兴趣区域格式
**典型数据：** 几何 ROI、选择
**用例：** ImageJ/Fiji 分析工作流程
**Python 库：**
- `read-roi`：`read_roi.read_roi_file('file.roi')`
- `roifile`：ROI 操作
**EDA 方法：**
- ROI类型分析（矩形、多边形等）
- 坐标提取
- ROI 属性（面积、周长）
- 群体分析（投资回报率集）
- Z 位置和时间信息

### .zip（ROI 集）
**描述：** ImageJ ROI 的 ZIP 存档
**典型数据：** 多个 ROI 文件
**用例：** 批量 ROI 分析
**Python 库：**
- `read-roi`：`read_roi.read_roi_zip('file.zip')`
- 标准`zipfile`模块
**EDA 方法：**
- 集合中的 ROI 计数
- ROI类型分布
- 空间分布
- 重叠 ROI 检测
- 命名约定

### .ome.tif / .ome.tiff - OME-TIFF
**描述：** 带有 OME-XML 元数据的 TIFF
**典型数据：** 具有丰富元数据的标准化显微镜
**用例：** 生物格式兼容存储
**Python 库：**
- `tifffile`：OME-TIFF 支持
- `AICSImageIO`：OME 读取
- `python-bioformats`：生物格式集成
**EDA 方法：**
- OME-XML 验证
- 物理尺寸提取
- 通道命名和波长
- 平面位置（Z、C、T）
- 仪器元数据
- 生物格式兼容性

### .ome.zarr - OME-ZARR
**描述：** ZARR 上的 OME-NGFF 规范
**典型数据：** 用于生物成像的下一代文件格式
**用例：** 云原生成像、大型数据集
**Python 库：**
- `ome-zarr-py`：官方实施
- `zarr`：底层数组存储
**EDA 方法：**
- 多尺度分辨率级别
- 元数据符合 OME-NGFF 规范
- 坐标变换
- 标签和 ROI 处理
- 云存储优化
- 块访问模式

### .klb - 凯勒实验室区块
**描述：** 用于大数据的快速显微镜格式
**典型数据：** 光片显微镜，延时
**用例：** 高通量成像
**Python 库：**
- `pyklb`：KLB 读写
**EDA 方法：**
- 压缩效率
- 块结构
- 多分辨率支持
- 阅读性能基准测试
- 元数据提取

### .vsi - 整个幻灯片成像
**描述：** 虚拟幻灯片格式（多个供应商）
**典型数据：** 病理切片，大马赛克
**用例：** 数字病理学
**Python 库：**
- `openslide-python`：多格式 WSI
- `tiffslide`：纯Python替代方案
**EDA 方法：**
- 金字塔等级计数
- 下采样因子
- 相关图像（宏、标签）
- 瓷砖尺寸和重叠
- MPP（微米每像素）
- 背景检测
- 组织分割

### .ndpi - 滨松 NanoZoomer
**描述：**滨松幻灯片扫描仪格式
**典型数据：** 全玻片病理图像
**用例：** 数字病理工作流程
**Python 库：**
- `openslide-python`：NDPI 支持
**EDA 方法：**
- 多分辨率金字塔
- 镜头和客观信息
- 扫描区域和放大倍数
- 焦平面信息
- 组织检测

### .svs - Aperio ScanScope
**说明：** Aperio 整个幻灯片格式
**典型数据：** 数字病理切片
**用例：** 病理图像分析
**Python 库：**
- `openslide-python`：SVS 支持
**EDA 方法：**
- 金字塔结构
- MPP校准
- 标签和宏观图像
- 压缩质量
- 缩略图生成

### .scn - 徕卡 SCN
**描述：** 徕卡幻灯片扫描仪格式
**典型数据：** 全玻片成像
**用例：** 数字病理学
**Python 库：**
- `openslide-python`：SCN 支持
**EDA 方法：**
- 瓦片结构分析
- 收集组织
- 元数据提取
- 放大倍数