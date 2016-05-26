# BVAR



### Рабочая цепочка:


1. Скачиваем данные

2. Удаляем сезонность в Eviews + ручное копирование в Excel 

На выходе:

- `/data/adjusted_data_x.xlsx`

3. Устанавливаем необходимые пакеты R, запускаем `/R/100_install_packages.R`

4. Делаем перевод данных в R, запускаем `/R/200_load_after_eviews.R`.

На входе: 

- `/data/adjusted_data_x.xlsx`

На выходе:

- `/data/df_2015_final.csv`: отфильтованы данные начиная с 1995 года, взяты логарифмы
    
- `/data/data_2015_after_eviews.csv`: просто конвертация исходного `xlsx`-файла
  
5. Алгоритм Банбуры для разных наборов переменных, запускаем `/R/550_table_replicator.R`

Работает долго.

Для работы нужны:

- `/R/500_replicate_banbura_function.R`
- `/R/400_model_funs.R`
- `/R/400_model_lists.R`
- `/R/500_banbura_funs.R`

На входе:

- `/data/df_2015_final.csv`

На выходе:

- файлы `.... forecasts_(variable name).Rds`: прогнозы по каждой переменной
- файлы `.... model_info_(variable name).Rds`: информация о моделях по каждой переменной
- `tables_rmsfe_raw.Rds`: данные для построения таблиц из работы

Todo: сделать автоматизацию для воспроизведения обеих таблиц без ручной корректировки

6. Получение широких таблиц из работы, запускаем `/R/550_raw_table_transform.R`

На входе:

- `/estimation/tables_rmsfe/tables_rmsfe_raw.Rds`: данные для построения таблиц из работы

На выходе:

- `/estimation/tables_rmsfe/table_rmsfe_final.csv`: табличка из работы

Todo: сделать автоматизацию для воспроизведения обеих таблиц без ручной корректировки


7. Нахождение оптимального множества моделей, model selection set, запускаем `/R/600_model_selection_set.R'

Работает долго.

На входе:

- файлы `.... forecasts_(variable name).Rds`: прогнозы по каждой переменной
- файлы `.... model_info_(variable name).Rds`: информация о моделях по каждой переменной

На выходе:

- файлы `.... loss_(variable name).Rds`: значение функции потерь для каждой переменной
- файлы `.... best_models_(variable name).Rds`: множество наилучших моделей для каждой переменной

Todo: при каком варианте гиперпараметров делаем (или с двумя как в случае таблиц?)


8. Описание множеств наилучших моделей, запускаем `/R/660_describe_mcs.R`

Графики и таблички для описания наилучших множеств.

На входе:

- файлы `.... best_models_(variable name).Rds`: множество наилучших моделей для каждой переменной

Todo: при каком варианте гиперпараметров делаем (или с двумя как в случае таблиц?)


### Тестирование 

#### Тестирование функции `replicate_banbura()`

1. Запускаем `550_test_replicate_banbura_function.R`

Для работы нужны:

- `/R/500_replicate_banbura_function.R`
- `/R/400_model_funs.R`
- `/R/400_model_lists.R`
- `/R/500_banbura_funs.R`

На входе:

- `/data/df_2015_final.csv`

### Разное

#### Попытка удалять сезонность в R

1. Переводим данные из Excel в R: `/R/200_load.R`

На входе: 

- `data/data_bvar_2015.xlsx`

На выходе:

- `/data/data_2015.csv`

2. Удаляем сезонность в R `/R/300_clean.R`:

На входе:

- `/data/data_2015.csv`

На выходе:

- `/data/df_2015_final.csv`

- `/data/df_2015_sa.csv`


#### Для загрузки ряда gas на www.seasonal.website

1. Отбираем нужное, `R/301_gas_problem.R`

На входе:

- `/data/data_2015.csv`

На выходе:

- `/data/gas_only.csv`

#### Попытка сравнения с carriero (не окончено)

1. Оцениваем нашим способом. Похоже параметры не те, `/R/350_matlab_carriero_test.R`

На входе:

- `/data/usa_data.mat`

### Эксперименты с model selection set

1. Первые шаги с model selection set, `/R/600_mcs.R`. 

2. Параллельные вычисления для model selection set, `/R/610_mcs_parallel.R`

3. Зависит ли model selection set от третьих альтернатив и как у него с ростом врени, `/R/610_mcs_3_alternative.R`

Увы: зависит от третьих альтернатив. Увы: время растёт квадратично по числу моделей.


