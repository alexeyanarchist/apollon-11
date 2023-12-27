# apollon 11

Используемые библиотеки
В данном проекте используются библиотеки math (для математических функций необходимых в расчетах, таких как sin и cos), prettytable (для построения красивой текстовой таблицы на основе полученных данных), matplotlib (для построения графиков на основе расчитанных данных) и numpy для корректной работы matplotlib.

Команда для установки всего необходимого:
pip install prettytable numpy matplotlib

Структура проекта
На данный момент проект разбит на 3 файла:

calculations.py - часть программы отвечающая за расчеты и получение таблицы.

graphics.py - часть программы, отвечающая за графическое представление расчитанных данных

data.py - часть программы, содержащая в себе все константы для расчетов

Инструкция к запуску
При запуске calc.py будут созданы или перезаписаны файлы table.csv и table.txt, а также будет выведена таблица, аналогичная той, которая находится в table.txt.

Примечание: файл table.csv из-за особенносьей кодировки некорректно открывается в некоторых табличных редакторвх, для корректной работы рекомендуется использовать excel

При запуске graphics.py будут созданы графики, графики будут сохранены в локальную папку graphic_img/ (в случае отсутсвия таковой она будет создана автоматически), помимо этого графики будут выведены через GUI
