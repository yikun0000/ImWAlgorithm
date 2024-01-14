#' 形态学重建
#'@description
#'
#' 这个函数执行形态学重建，使用指定的标记和掩膜图像，以及给定的结构元素。
#'
#' @param marker 标记图像，通常是腐蚀或变换后的图像。
#' @param mask 掩膜图像，定义了重建的极限。
#' @param kernel2 结构元素的大小，用于定义膨胀操作的范围。
#'
#' @return 返回形态学重建后的图像。
#' @name morphoReconstruct 
#' @examples 
#' # 假设有两个图像 marker 和 mask
#' # 运行形态学重建
#' result <- morphoReconstruct(marker, mask, 5)
#' @export 
morphoReconstruct <- function(marker, mask,kernel2) {
  old <- marker
  repeat {
    new <- dilate(old, makeBrush(kernel2, shape = 'diamond')) & mask
    if (all(new == old)) {
      break
    }
    old <- new
  }
  return(new)
}

#' 处理图像
#' @description
#' 
#' 这个函数对输入的彩色图像执行一系列的形态学操作，包括梯度计算、闭合、重建和分水岭分割。
#'
#' @param image 输入的彩色图像。
#' @param kernel1 计算形态学梯度的核大小。
#' @param kernel2 形态学闭重构的核大小。
#' @param h 执行H-minima变换的高度阈值。
#' @param d0 二阶巴特沃斯低通滤波器的截止频率。
#'
#' @return 返回处理后的图像列表，包括旋转后的图像和分水岭结果。
#' @name processImage
#' @examples 
#' # 假设 image 是一个加载的彩色图像
#' # 处理图像并获取结果
#' result <- processImage(image, 5, 5, 1, 30)
#' # 查看分水岭结果
#' plot(result$watershed_result)
#' @export 
processImage <- function(image,kernel1 = 5,kernel2 = 5, h = 1,d0 = 30) {
  #image为输入的彩色图像,kernel1为计算形态学梯度的核大小,kernel2为形态学闭重构的核大小
  #h 为执行H-minima的高度
  #order 为二阶巴特沃斯低通滤波的参数
  # 应用高斯模糊
  filtered_image <- gblur(image, sigma=1)
  
  # 执行形态学膨胀
  dilated <- dilate(filtered_image, makeBrush(kernel1, shape = 'diamond'))
  
  # 获取图像的尺寸
  row1 <- dim(image)[1]
  col1 <- dim(image)[2]
  
  # 执行形态学腐蚀
  eroded <- erode(filtered_image, makeBrush(kernel1, shape = 'diamond'))
  
  # 计算形态学梯度
  morphological_gradient <- dilated - eroded
  
  # 转置并获取最大梯度
  R <- t(matrix(morphological_gradient[,,1], nrow = row1, ncol = col1))
  G <- t(matrix(morphological_gradient[,,2], nrow = row1, ncol = col1))
  B <- t(matrix(morphological_gradient[,,3], nrow = row1, ncol = col1))
  max_gradient <- pmax(pmax(R, G), B)
  
  # 应用巴特沃斯低通滤波
  filtered_image <- butterworthLowPassFilter(max_gradient,d0, 2)
  filtered_image[filtered_image < 0] <- 0
  
  # 创建结构元素
  structuring_element <- makeBrush(kernel2, shape = 'diamond')
  
  # 执行形态学闭操作
  closed_image <- closing(filtered_image, structuring_element)
  
  
  # 应用形态学腐蚀
  filtered_image1 <- round(closed_image * 255)
  
  eroded <- erode(filtered_image1, makeBrush(5, shape = 'diamond'))
  
  # 应用形态学重建
  reconstructed <- morphoReconstruct(eroded, filtered_image1,kernel2)
  
  # 计算重建开操作
  opened <- opening(filtered_image, makeBrush(kernel2, shape = 'diamond'))
  
  # 找到极小值区域
  minima <- (filtered_image1 - reconstructed) >= h
  
  
  para1 = pmax(minima*max(filtered_image1),filtered_image1)
  # 应用分水岭算法
  watershed_result <- watershed(para1)
  
  # 旋转分水岭结果
  rotate_image <- rotate(watershed_result, angle = 90)
  flopped_image <- flop(rotate_image)
  # 返回处理后的图像，可以根据需要返回其他中间结果或最终结果
  return(list(rotate_image = flopped_image, watershed_result = watershed_result))
}














