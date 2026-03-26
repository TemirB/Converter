# Makefile для компиляции программы

# Путь к установке ROOT
ROOT_DIR = /home/tahea/root_install

# Путь к библиотекам McDst
MCDST_DIR = /home/tahea/external/McDst/build
MCDST_LIB = /home/tahea/external/McDst/include

# Компилятор
CXX = g++

# Флаги компилятора
CXXFLAGS = $(shell root-config --cflags) -fPIC -Wall -std=c++17
CXXFLAGS += -I$(MCDST_LIB) -I$(ROOT_DIR)/include -Iinclude
LIBS = $(shell root-config --glibs) -L$(ROOT_DIR)/lib -lCore -lRIO -lTree -lHist -lGraf -lGpad -lPhysics -lThread -lm
LIBS += -L$(MCDST_DIR) -lMcDst

# Имя исполняемого файла
TARGET = converter

# Директории
SRC_DIR = src
BUILD_DIR = build
BUILD_DIR_PROFILE = build_profile

# Источники
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
OBJS_PROFILE = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR_PROFILE)/%.o, $(SRCS))

# Флаги для профайлинга (gprof)
PROFILE_FLAGS = -pg

# Создать директории при необходимости
$(shell mkdir -p $(BUILD_DIR))
$(shell mkdir -p $(BUILD_DIR_PROFILE))

# Правила сборки
all: $(TARGET)

# Сборка исполняемого файла (обычная)
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

# Сборка с профайлингом (gprof)
profile: $(BUILD_DIR_PROFILE)/$(TARGET)
	cp $(BUILD_DIR_PROFILE)/$(TARGET) $(TARGET)_profile

# Профилированная версия исполняемого файла
$(BUILD_DIR_PROFILE)/$(TARGET): $(OBJS_PROFILE)
	$(CXX) -pg -o $@ $^ $(LIBS)

# Компиляция объектных файлов (обычная)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Компиляция объектных файлов (с профайлингом)
$(BUILD_DIR_PROFILE)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(PROFILE_FLAGS) -c $< -o $@

# Запуск профайлера (автоматический анализ)
run-profile: profile
	./$(TARGET)_profile
	gprof $(TARGET)_profile gmon.out > profile_report.txt
	@echo "Профайлинг завершён. Отчёт в profile_report.txt"

# Очистка
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Полная очистка (включая профили)
distclean: clean
	rm -rf $(BUILD_DIR_PROFILE) $(TARGET)_profile gmon.out profile_report.txt
	rm -f *~

# Помощь
help:
	@echo "Доступные цели:"
	@echo "  make          - Обычная сборка"
	@echo "  make profile  - Сборка с профайлингом (создаёт converter_profile)"
	@echo "  make run-profile - Сборка + запуск + анализ gprof"
	@echo "  make clean    - Очистка обычных файлов"
	@echo "  make distclean - Полная очистка (включая профиль)"