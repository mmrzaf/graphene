CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinclude
AR = ar
ARFLAGS = rcs

SRC_DIR = src
INC_DIR = include
TEST_DIR = tests
BUILD_DIR = build
LIB_DIR = lib

SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SOURCES))
LIBRARY = $(LIB_DIR)/libgraphene.a

TEST_SOURCES = $(wildcard $(TEST_DIR)/*.c)
TEST_BINARIES = $(patsubst $(TEST_DIR)/%.c,$(BUILD_DIR)/%,$(TEST_SOURCES))

.PHONY: all clean test dirs

all: dirs $(LIBRARY)

dirs:
	@mkdir -p $(BUILD_DIR) $(LIB_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBRARY): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^
	@echo "Built library: $(LIBRARY)"

$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.c $(LIBRARY)
	$(CC) $(CFLAGS) $< -L$(LIB_DIR) -lgraphene -o $@

test: dirs $(LIBRARY) $(TEST_BINARIES)
	@echo "Running tests..."
	@for test in $(TEST_BINARIES); do \
		echo ""; \
		./$$test || exit 1; \
	done
	@echo ""
	@echo "All tests passed!"

clean:
	rm -rf $(BUILD_DIR) $(LIB_DIR)
	rm -f test_temp.edgelist

install: $(LIBRARY)
	@mkdir -p /usr/local/lib /usr/local/include/graphene
	cp $(LIBRARY) /usr/local/lib/
	cp $(INC_DIR)/*.h /usr/local/include/graphene/
	@echo "Installed to /usr/local"

