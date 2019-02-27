本项目是由cmake构建的cpp，项目中包含简略的cmake构建文件


依赖：opencv library和ffmpeg，需要预先装好相关库
默认相关头文件和库在cmake中：
set(INC_DIR /usr/local/include)
set(LINK_DIR /usr/local/lib)
需要根据当前文件夹进行修改

确定完成后可以用cmake进行项目构建，用法：
./FrameExtractor -i in_video -o out_dir(默认./out_put) -S(开启按秒随机抽帧)