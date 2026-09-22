CC     ?= gcc
CFLAGS ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Iinclude
BUILD  := build

.PHONY: all test sim clean

all: test sim

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/foctest: src/foc_math.c src/foc_ctrl.c test/test_foc.c include/foc.h | $(BUILD)
	$(CC) $(CFLAGS) src/foc_math.c src/foc_ctrl.c test/test_foc.c -o $@ -lm

$(BUILD)/focsim: src/foc_math.c src/foc_ctrl.c src/main.c include/foc.h | $(BUILD)
	$(CC) $(CFLAGS) src/foc_math.c src/foc_ctrl.c src/main.c -o $@ -lm

test: $(BUILD)/foctest
	./$(BUILD)/foctest

sim: $(BUILD)/focsim
	./$(BUILD)/focsim pmsm
	./$(BUILD)/focsim bldc
	./$(BUILD)/focsim ff
	./$(BUILD)/focsim fw

clean:
	rm -rf $(BUILD)
