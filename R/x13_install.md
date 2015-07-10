Установка X13 под mac

1. Скачать исходники файл `x13assrc.tar.gz` c [census.gov](http://www.census.gov/srd/www/x13as/x13down_unix.html)
2. Разархивировать. Получится папка с кучей файлов на фортране.
3. Запустить терминал
4. Набрать `which gfortran`
скопировать полученный адрес нахождения фортрана (у меня получилось `/usr/local/bin/gfortran`)
5. отредактировать файл `makefile.g77`
В первых двух строчках исправить путь на полученный в пункте 4.
например:
FC        = /usr/local/bin/gfortran
LINKER    = /usr/local/bin/gfortran
6. В терминале набрать `make -f makefile.g77`
7. В R:
```r
Sys.setenv(X13_PATH ="путь к папке с файлами")
library("seasonal")
checkX13()
```
