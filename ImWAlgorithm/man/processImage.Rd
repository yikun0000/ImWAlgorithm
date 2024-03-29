% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Watershed.R
\name{processImage}
\alias{processImage}
\title{处理图像}
\arguments{
\item{image}{输入的彩色图像。}

\item{kernel1}{计算形态学梯度的核大小。}

\item{kernel2}{形态学闭重构的核大小。}

\item{h}{执行H-minima变换的高度阈值。}

\item{d0}{二阶巴特沃斯低通滤波器的截止频率。}
}
\value{
返回处理后的图像列表，包括旋转后的图像和分水岭结果。
}
\description{
这个函数对输入的彩色图像执行一系列的形态学操作，包括梯度计算、闭合、重建和分水岭分割。
}
\examples{
# 假设 image 是一个加载的彩色图像
# 处理图像并获取结果
result <- processImage(image, 5, 5, 1, 30)
# 查看分水岭结果
plot(result$watershed_result)
}
