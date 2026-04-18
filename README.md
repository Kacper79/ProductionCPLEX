# ProductionCPLEX
To repozytorium implementuje w narzedziu CPLEX (w Concert API z wykorzystaniem C++) problem produkcji (Closed-Loop Supply Chain) z [*artykulu*](https://link.springer.com/article/10.1007/BF03342742).
Artykul opisuje bazowy model CLSP-RM-SS realizujący problem - multi-product lot-sizing problem with product returns and remanufacturing subject to a capacity constraint oraz przedstawia różne warianty modeli produkcyjnych.
W tym repozotorium jest zrealizowany model CLSP-RM-SS (bazowy) oraz problem CLSP-RM (problem górny SPP + problem dolny SLULSP-RM) z wykorzystaniem techniki generacji kolumn do generowania poprawnego wzorca (schedules). 

# Pobieranie
```git clone https://github.com/Kacper79/DystroCPLEX.git```

```cd DystroCPLEX```

# Budowanie
Wazne jest aby w pliku CMakeLists.txt ustawic lokalnie prawidlowa sciezke do folderu gdzie zostal zainstalowany CPLEX na komputerze.

```set(CPLEX_ROOT_PATH "nowa_sciezka_do_folderu_cplex") ```
gdzie sciezka (na Windowsie) musi byc podawana ze slashami zwyklymi (/) a nie odwrotnymi (\\) ktore sa problematyczne dla CMake

Jezeli wersje CPLEX sie roznia, trzeba bedzie takze zmienic nazwe biblioteki w pliku CMakeLists przy podpinaniu bibliotek (.dll/.so) np. zmienic na starsza wersje z cplex2212 na cplex2211 (jak w pracowni UKSW)

Potem wystarczy zbudowac za pomoca CMake
```cmake -S . -B build```

```cmake --build build --config Release```
