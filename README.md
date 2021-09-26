# Структура репозитория
### methods_approximation_monitoring_data
* **[Classic_LS_SW](https://github.com/aleksandra-shchegoleva/matlab/blob/master/methods_approximation_monitoring_data/Classic_LS_SW.m)** - подбор коэффициентов для дискретной системы "хищник-жертва" с помощью классического метода наименьших квадратов
* **[SW_2019_data](https://github.com/aleksandra-shchegoleva/matlab/blob/master/methods_approximation_monitoring_data/SW_2019_data)** - шведские данные мониторинга за 2019 год (некоторые виды)
### NAS
* **[plot_NAS](https://github.com/aleksandra-shchegoleva/matlab/blob/master/NAS/plot_NAS.m)** - код для построения системы с применением алгоритма NAS
* **[apply_NAS](https://github.com/aleksandra-shchegoleva/matlab/blob/master/NAS/apply_NAS.pdf)** - файл с описанием системы
### Файлы, лежащие вне папок
Здесь лежат общие файлы без привязки к каким-либо данным
* **[System_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/System_U.m)** - код для построения системы с управлением
* **[Stability_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Stability_U.m)** - код, находящие коэффициенты системы с управлением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **`System_B`** - код, строющий график системы с управлением и ограничением
* **`Stability_B`** - код, находящие коэффициенты системы с управлением и ограничением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **`System_B_r`** - код, строющий график системы с управлением, ограничением и коэффициентом r
* **`Stability_B_r`** - код, находящие коэффициенты системы с управлением, ограничением и коэффициентом r, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **[System_U_Runge_Kutta_4](https://github.com/aleksandra-shchegoleva/matlab/blob/master/System_U_Runge_Kutta_4.m)** - численное решение системы "хищник-жертва с питанием" методом Рунге-Кутта 4 порядка
### ModelApp
Здесь лежат файлы для программы построения графиков для системы с управлением с графическим интерфейсом.

---
Тестовые значения увеличения различных данных
| Станция        | Популяция           | Биомасса  |
| ------------- |:-------------:| -----:|
| 7      | 5 | 0.5 |
| 2      | 1      |   0.5 |
| 12 | 2      |    0.5 |


* **[parse_excel](https://github.com/aleksandra-shchegoleva/matlab/blob/master/ModelApp/parse_excel.m)** - код, читающий данные о популяции и биомассе из Excel файла (лежит в этой же папке)
* **[parse_date](https://github.com/aleksandra-shchegoleva/matlab/blob/master/ModelApp/parse_date.m)** - код, читающий данные о времени и станциях из Excel файла (лежит в этой же папке)
* **[plot_system_full](https://github.com/aleksandra-shchegoleva/matlab/blob/master/ModelApp/plot_system_full.m)** - построение системы с управлением по введенным начальным данными
* **[app.mlapp](https://github.com/aleksandra-shchegoleva/matlab/blob/master/ModelApp/app.mlapp)** - файл App Designer
### Real_data
Здесь лежат файлы, использующие данные мониторинга из файла Excel ([формат файла](https://drive.google.com/file/d/1T4Fsw-0qFkj_fwCRzScHDM4GHsqTIcmT/view?usp=sharing))
* **[System_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_U.m)** - код, строющий график системы с управлением, исходные данные берутся из файла Excel
* **[Search_coeff](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/Search_coeff.m)** - код, находящие коэффициенты для системы без управления, при которых система наиболее приближена к данным мониторинга ([описание способа](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/Method.txt))
* **[Stability_U](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/Stability_U.m)** - код, находящие коэффициенты системы с управлением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
* **[plot_system_full](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/plot_system_full.m)** - построение системы с управлением по введенным начальным данными (объединение файлов Search_coeff, Stability_U, System_U)
* **[System_B_plot](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_B_plot.m)** - код, строющий график системы с управлением и ограничением
* **[System_B](https://github.com/aleksandra-shchegoleva/matlab/blob/master/Real_data/System_B.m)** - код, находящие коэффициенты системы с управлением и ограничением, при которых она стабильная и достигается цель, на выходе получаем массив значений коэффициентов
