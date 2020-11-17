# Структура репозитория
### Файлы, лежащие вне папок
Здесь лежат общие файлы без привязки к каким-либо данным
* **[System_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/System_U.m)** - код, строющий график системы с управлением
* **[Stability_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Stability_U.m)** - код, находящие коэффициенты системы с управлением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **`System_B`** - код, строющий график системы с управлением и ограничением
* **`Stability_B`** - код, находящие коэффициенты системы с управлением и ограничением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **`System_B_r`** - код, строющий график системы с управлением, ограничением и коэффициентом r
* **`Stability_B_r`** - код, находящие коэффициенты системы с управлением, ограничением и коэффициентом r, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
### Real_data
Здесь лежат файлы, использующие данные мониторинга из файла Excel ([формат файла](https://drive.google.com/file/d/1T4Fsw-0qFkj_fwCRzScHDM4GHsqTIcmT/view?usp=sharing))
* **[System_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_U.m)** - код, строющий график системы с управлением, исходные данные берутся из файла Excel
* **[Search_coeff](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/Search_coeff.m)** - код, находящие коэффициенты для системы без управления, при которых система наиболее приближена к данным мониторинга (описание способа)
* **[Stability_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/Stability_U.m)** - код, находящие коэффициенты системы с управлением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **[plot_system_full](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/plot_system_full.m)** - построене системы с управлением по введенным начальными данными (объединение файлов Search_coeff, Stability_U, System_U)
* **[System_B_plot](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_B_plot.m)** - код, строющий график системы с управлением и ограничением
* **[System_B](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_B.m)** - код, находящие коэффициенты системы с управлением и ограничением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
