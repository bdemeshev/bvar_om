bmr
* Normal-inverse-Wishart Prior
* Minnesota Prior
  * (Koop and Korobilis, 2010, p. 7) version
  * (Canova, 2007, p. 378) version
    * harmonic decay d(l)=l^{H_4}
    * geometric decay d(l)=H_4^{1-l}
* steady state prior


Eviews:
* Litterman/Minnesota prior. (koop korobilis version) v
* Normal-Wishart prior. v
* Sims-Zha Normal-Wishart prior.
* Sims-Zha Normal-flat. v
* Diffuse prior 
? conjugate Normal-Wishart 


http://www.ihs.ac.at/vienna/resources/Economics/Papers/20120426_Slides_Giannone.pdf
* Minnesota prior (Litterman, 1980 and 1986) 
* Inverse-Wishart prior on \Sigma
* Sum-of-coefficients prior (Doan, Litterman and Sims 1984) 
* Single-unit-root prior (Sims, 1993)

http://personal.strath.ac.uk/gary.koop/kk3.pdf
Gary Koop (from bvar_full)
* Diffuse (Jeffreys) (MC integration)
* Minnesota (MC integration)
* Normal-Wishart (MC integration)
* Independent Normal-Wishart (Gibbs)
* SSVS in mean - Wishart (Gibbs)
* SSVS in mean - SSVS in in covariance (Gibbs)
SSVS -- Stochastic Search Variable Selection

y_t=a_0+A_1 y_{t-1} + \varepsilon_t


http://www.ecb.europa.eu/events/pdf/conferences/schorfheide.pdf
* Priors from General Equilibrium Models for VARs

Обозначения:
y_t=\Phi_0+\Phi_1 y_{t-1} +u_t



msbvar:

szbvar:
* Sims-Zha Normal-Wishart prior
* Sims-Zha Normal-flat prior
* (?) flat-flat prior

szbsvar
* "model described by Sims and Zha (1998) and Waggoner and Zha (2003)"


уже начатая статья
y_t=v+A_1 y_{t-1} + \ldots + u_t





Sims, C.A. and Tao A. Zha. 1998. "Bayesian Methods for Dynamic Multivariate Models." International Economic Review. 39(4):949-968.

Waggoner, Daniel F. and Tao A. Zha. 2003a. "A Gibbs sampler for structural vector autoregressions" Journal of Economic Dynamics \& Control. 28:349–366.

Waggoner, Daniel F. and Tao A. Zha. 2003b. "Likelihood preserving normalization in multiple equation models". Journal of Econometrics. 114: 329–347.

Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian Time Series Modeling and the Study of Politics: Theory Testing, Forecasting, and Policy Analysis". Political Analysis 14(1):1-36.


Протокол:

1.       Сезонная корректировка
2.       Взятие логарифмов для всех переменных, исключая те, что выражены  в процентах.
3.       Проверка на стационарность для определения априорного  матожидания (опционально. Carriero этого вообще не делает) 
4. Формируем выборки с 3, 5, 7, 15 и 23 переменными 
5. Оценивать нужно RW, VAR(3)=BVAR(\lambda=\infty), VAR(5), BVAR_ Minnesota, сопр. N-IW
6.  Количество лагов для  BVAR выбираем путем максимизации marginal density. Лаги смотрим от 1 до 12. Первая точка в обучающей выборке – январь 1996, количество исходных наблюдений – 150, выборка не меняется с изменением числа лагов. Второй вариант – количество исходных наблюдений – 200. Из любопытства можно сравнить выбранное предпочитаемое количество лагов с результатом минимизации информационных критериев.  
7.  Подсчет \lambda по банбуровской схеме.  Процентную ставку во внимание не принимаем (хотя, как мы выяснили, она не вредная, но и не полезная)
8 Подсчет \lambda_1 путем максимизации marginal density. Можно по отдельности, можно максимизировать одновременно по \lambda_1 и по p.  \lambda_2  надо ставить равным 1. А по \lambda_3 и, возможно,  \lambda_4 можно оптимизировать. Вопрос, хватает ли для этого данных? 
9. Прогнозы строим на 1,3, 6, 9 и 12 месяцев. Считаем темп роста прогнозируемых переменных (можно начать с прироста IP   и инфляции, потом, если результаты будут приятными, можно расширить число прогнозируемых переменных). 
10. Сдвигаем наблюдения на единицу и делаем то же, что в п. 9. Количество наблюдений в выборке всегда 150. Последнее прогнозируемое наблюдение – апрель 2015, количество прогнозов на 1 период будет больше, чем на 3, а прогнозов на 3 будет больше, чем на 6
11. Считаем MSFE по всем прогнозам, делаем выводы. 
12.  Потом решаем, что были глубоко неправы везде, где можно. Возвращаемся и переделываем. 


о матрицах:
- svd, Schur, QR, ...
- разобраться с Муром-Пенроузом псевдообратной (решить пару примеров, свойства)
- разобраться с разложением Холецкого (решить пару примеров, свойста)

о классификации:
- отличие Сим-Жа от сопряженного нормального-Вишарта? похоже, что Сим-Жа более общий случай
- при каких априорных параметрах апостериорное среднее совпадёт с ML VAR?
будет ли это неинформативное-Джеффри?


мелочи:
- у банбуры пока идет речь о подборе лямбды все msfe --- in-sample

дизайн пакета:
priors <- bvar_priors(type="conjugate", style/spec="carriero", Y, p=4, lambdas=c(.....))
model <- bvar_estimate(priors)
summary(model)
пример тестов из coda:

forecast(model, h=)
forecast(model, out-of-sample=FALSE)


потестить:
- bvar optimal и var-3 (ок для fast_forecast)
- без fast_forecast
- при добавлении искуственных наблюдений (просто, а также sc, io) какое T брать при 
преобразовании v_post=v_prior + T

19 августа 2015
18:51 возможная причина ошибки: у BVAR-3 n_lag=1 omsfe резко растут при переходе с h=1 до h=2 (ok, поймали ошибку: не надо копировать прошлый x_t при h>1)
19:17 добавить названия строк/столбцов в colnames/rownames Phi (ok)
сделать colnames/rownames Omega

20 августа 2015
13:01 сделать оценку VAR с автоподбором лагов (ok)
13:02 сделать табличку как у банбуры стр. 79 (ok)
18:09 сделать так, чтобы пользователь мог не вносить в априорный список вспомогательные 
элементы
http://stackoverflow.com/questions/7719741/how-to-test-if-list-element-exists
18:45 сверить теоретическое апостериорное среднее и выборочное --- (ок)

21 августа 2015
0:25 сравнить средние коэффицинты: аналитическое апостериорное, выборочное среднее при двух способах генерировать коэффициенты (ok), sample sd (попалась разница)
затем сравнить msfe
сделать код для многих ядер процессора: в самом bvar_conjugate и в estimate
почистить код на минимальный priors (чтобы минимальный работал)
16:02 у способов svd/chol отличаются sample sd of Phi post (!corrected,   не там было транспонирование)
16:18  может 13 лагов лучше чем 12?

23 августа 2015
11:26 изменить выборку (убрать 1995 год, например)
11:27 зафиксировать лямбды на base scenario Carriero
11:27 три части выборки: оценка, подбор лямбды под оптим. прогнозы, оценка
11:31 ввести id процесса оценивания, затем читалку результатов оценивания
14:30 гармония явного подхода с заданием априорных матриц: dummy_sс + dummy_io влияют только на априорные матрицы!!! но не включаются более никуда!!!!
14:55 больше наборов данных для сравнения прямо в пакет

24 августа 2015
* Обнаружена ошибка в forecast rw in-sample, бралось время t=1, а надо брать с момента предсказания минус один. Исправлено.
* Еще одна ошибка в forecast rw in-sample. Делался прогноз от первой точки на h шагов вперед, а надо было от предыдущей на один шаг вперед. Исправлено.
* задача на сегодня: разобраться с гармонией явного подхода с заданием априорных матриц.
* fit_set в bvar_out_list (есть, ок)

25 авгутса 2015
* best_lambda наличие повторных строк (нет, ок)
* дизайн реализации двух подходов (ситуация слабо осложнена наличием nu всегда)
obs2prior
prior2obs (??? как выглядит решения для произвольных prior)
lambda2obs
lambda2hyper (двойная конвертация labmda -> obs -> prior)
estimate from obs (выбрать более робастную из from obs/prior)
estimate from prior
* в пакете сделать prior$Omega вместо prior$Omega_prior?
* хранить Omega или Omega^{-1}? 
* разделить симуляции по banbura от отчёта по банбуре? 
в отчете сверху указывается желаемое число лагов и прочее
* сделать корректные типы, чтобы убрать warnings() (ок)
* отказаться от длинных списков - оставить только широкие (ок) 

26 августа 2015
* стандартизировать test/verbose, добавить печатаемую инфу, добавить проверку состоятельности
* разделить оценку и файл с параметрами и отчетом ?
* добавить удобную оценивалку одной модели
* нужен ли constant в setup модели? при прогнозе вроде как да (да)
* нужен ли Y_in в setup? если он легко рековерится? (нет)

8 сентября 2015
* корректность mdd (marginal data density)? Сравнение с Карьеро, рассмотрение предельных случаев, где можно руками посчитать
* неединственность оптимальной mdd (ок)
* бесконечные лямбда --- протестировать все функции
* дописать rmsfe (ok)
* сделать разные v_prior в списках моделей 

9 сентября 2015
* выбор первой попавшейся оптимальной mdd в случае равенства, 
стандартизировано с replicate banbura
* rmsfe for mdd ok
* в чем все (!) отличия кода каррьеро от нашего?

10 сентября 2015
* СТРАННАЯ ошибка: source() в R работает, пошаговое исполнение в Rstudio - ok
source() в Rstudio не работает
"Error in eval(expr, envir, enclos) : object at index 356 not a data.frame"
" incompatible type (data index: 414, column: 't', was collecting: numeric (dplyr::Collecter_Impl<14>),"
свежий RRO и Rstudio и пакеты после переустановки
  - вымарать foreach? (не помогло) попробовать классику?
  - перезагрузка? (кажется помогло?)
* автоматизация номеров столбцов в конце
* (ok) проверка корректности оценок при sigma^2 из AR(1) 
* удобная функция для одной модели
* с каким omsfe сравнивать? с RW всегда? с моделью выбранной по тесту? с AR() с нужным коэффициентом?
- сравниваем с RW всегда, а установка delta влияет только на prior
 

