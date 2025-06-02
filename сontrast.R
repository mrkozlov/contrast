# Установка и загрузка необходимых пакетов 
library(sf)
library(terra)
library(stats)
library(writexl)

# Очистка памяти
gc()

# Настройки terra для работы с диском
temp_dir <- "D:/temp_terra6"
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = temp_dir, memfrac = 0.3)

# Загрузка точек присутствия
# Предполагается, что train_balanced уже в памяти как sf-объект
# Если точки в shapefile: train_balanced <- st_read("train_balanced.shp")

# Установка CRS для train_balanced (WGS84)
st_crs(train_balanced) <- "+proj=longlat +datum=WGS84"
cat("CRS train_balanced:", st_crs(train_balanced)$input, "\n")

# Загрузка полигона
polygon <- st_read("15_11Layers.shp")

# Приведение CRS полигона к WGS84
if (st_crs(polygon) != st_crs(train_balanced)) {
  polygon <- st_transform(polygon, crs = "+proj=longlat +datum=WGS84")
  cat("CRS полигона приведено:", st_crs(polygon)$input, "\n")
}
polygon_vect <- vect(polygon)

# Загрузка растровых слоев
path_30s <- "t30"
rasters_30s <- list.files(path = path_30s, pattern = "\\.tif$", full.names = TRUE)

# Выбираем только указанные слои для маски SRE (точное совпадение)
selected_layers <- c("wc2.1_30s_bio_8", "wc2.1_30s_bio_14", "wc2.1_30s_bio_3")
raster_names <- gsub(".tif", "", basename(rasters_30s))
selected_indices <- which(raster_names %in% selected_layers)
if (length(selected_indices) != length(selected_layers)) {
  stop("Не удалось найти все указанные слои для маски. Проверьте имена файлов в папке t30.")
}
selected_rasters <- rasters_30s[selected_indices]
env_layers_mask <- rast(selected_rasters)

# Проверка количества выбранных слоев для маски
cat("Количество слоев для маски SRE:", nlyr(env_layers_mask), "\n")
cat("Имена слоев для маски:", names(env_layers_mask), "\n")

# Проверка и установка CRS растров для маски
if (is.na(crs(env_layers_mask)) || crs(env_layers_mask) != "+proj=longlat +datum=WGS84") {
  crs(env_layers_mask) <- "+proj=longlat +datum=WGS84"
  cat("CRS растров для маски установлено:", crs(env_layers_mask), "\n")
}

# Установка имен слоев для маски
names(env_layers_mask) <- gsub(".tif", "", basename(selected_rasters))

# Обрезка растров для маски по полигону
env_layers_mask <- crop(env_layers_mask, polygon_vect)

# Проверка и удаление существующего файла перед записью
cropped_file_mask <- file.path(temp_dir, "cropped_rasters_mask.tif")
if (file.exists(cropped_file_mask)) {
  result <- unlink(cropped_file_mask, force = TRUE)
  if (result == 0) {
    cat("Файл", cropped_file_mask, "успешно удален.\n")
  } else {
    stop("Не удалось удалить существующий файл", cropped_file_mask, ". Проверьте права доступа или закройте файл.")
  }
}
writeRaster(env_layers_mask, filename = cropped_file_mask, overwrite = TRUE)

# Проверка успешности записи
if (file.exists(cropped_file_mask)) {
  env_layers_mask <- rast(cropped_file_mask)
  gc()
  cat("Растры для маски успешно записаны в", cropped_file_mask, "\n")
} else {
  stop("Не удалось записать обрезанные растры для маски в файл")
}

# Извлечение значений экологических переменных для точек присутствия (для маски)
presence_vect <- vect(train_balanced)
env_presence_mask <- terra::extract(env_layers_mask, presence_vect)[, -1]
if (any(is.na(env_presence_mask))) {
  warning("NA обнаружены в данных. Удалено строк: ", sum(!complete.cases(env_presence_mask)))
  env_presence_mask <- na.omit(env_presence_mask)
  env_presence_mask <- env_presence_mask[valid_rows, ]
  presence_vect <- presence_vect[valid_rows]
}
gc()

# Вычисление квантилей (базовые 0.025)
quantiles <- apply(env_presence_mask, 2, function(x) {
  quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
})

# Оптимизированная функция для создания маски SRE
create_sre_mask <- function(env_layers, quantiles, temp_dir) {
  mask <- env_layers[[1]] * 0 + 1
  for (i in 1:nlyr(env_layers)) {
    layer <- env_layers[[i]]
    q_low <- quantiles[1, i]
    q_high <- quantiles[2, i]
    temp_mask <- layer < q_low | layer > q_high
    temp_file <- file.path(temp_dir, paste0("temp_mask_", i, ".tif"))
    if (file.exists(temp_file)) unlink(temp_file)
    writeRaster(temp_mask, filename = temp_file, overwrite = TRUE)
    mask <- mask * rast(temp_file)
    gc()
    cat("Частота значений в маске после слоя", i, ":\n")
    print(freq(mask))
  }
  return(mask)
}

# Создание маски SRE
sre_mask <- create_sre_mask(env_layers_mask, quantiles, temp_dir)

# Проверка маски перед записью
cat("Частота значений в маске SRE (до записи):\n")
print(freq(sre_mask))

# Сохранение маски на диск
sre_mask_file <- file.path(temp_dir, "sre_mask.tif")
if (file.exists(sre_mask_file)) unlink(sre_mask_file)
writeRaster(sre_mask, filename = sre_mask_file, overwrite = TRUE)
sre_mask <- rast(sre_mask_file)
gc()

# Проверка маски после записи
cat("Частота значений в маске SRE (после записи):\n")
print(freq(sre_mask))

# Замена 0 на NA после записи
sre_mask[sre_mask == 0] <- NA

# Проверка маски после замены
cat("Частота значений в маске SRE (после замены 0 на NA):\n")
print(freq(sre_mask))

# Ограничение маски полигоном
sre_mask <- mask(sre_mask, polygon_vect)

# Проверка маски после полигона
cat("Частота значений в маске после полигона:\n")
print(freq(sre_mask))

# Проверка, есть ли ненулевые пиксели
if (all(is.na(values(sre_mask)))) {
  stop("Маска SRE пуста (все значения NA). Попробуйте ослабить квантиль или проверить выбор слоев.")
}

# Генерация псевдоотсутствий (соотношение 1:1)
n_pseudo <- nrow(train_balanced)  # Соотношение 1:1
cat("Количество точек присутствия:", nrow(train_balanced), "\n")
cat("Количество псевдоотсутствий:", n_pseudo, "\n")
pseudo_absences <- spatSample(sre_mask, size = n_pseudo, method = "random", na.rm = TRUE, as.points = TRUE)

# Извлечение координат псевдоотсутствий
pseudo_coords <- crds(pseudo_absences)

# Загрузка всех 20 слоев для анализа
all_rasters <- rasters_30s
env_layers_all <- rast(all_rasters)

# Проверка количества всех слоев
cat("Количество всех слоев для анализа:", nlyr(env_layers_all), "\n")
cat("Имена всех слоев:", names(env_layers_all), "\n")

# Проверка и установка CRS для всех слоев
if (is.na(crs(env_layers_all)) || crs(env_layers_all) != "+proj=longlat +datum=WGS84") {
  crs(env_layers_all) <- "+proj=longlat +datum=WGS84"
  cat("CRS всех растров установлено:", crs(env_layers_all), "\n")
}

# Обрезка всех слоев по полигону
env_layers_all <- crop(env_layers_all, polygon_vect)

# Проверка и удаление существующего файла перед записью
cropped_file_all <- file.path(temp_dir, "cropped_rasters_all.tif")
if (file.exists(cropped_file_all)) {
  result <- unlink(cropped_file_all, force = TRUE)
  if (result == 0) {
    cat("Файл", cropped_file_all, "успешно удален.\n")
  } else {
    stop("Не удалось удалить существующий файл", cropped_file_all, ". Проверьте права доступа или закройте файл.")
  }
}
writeRaster(env_layers_all, filename = cropped_file_all, overwrite = TRUE)

# Проверка успешности записи
if (file.exists(cropped_file_all)) {
  env_layers_all <- rast(cropped_file_all)
  gc()
  cat("Все растры успешно записаны в", cropped_file_all, "\n")
} else {
  stop("Не удалось записать обрезанные растры для анализа в файл")
}

# Проверка типа и структуры pseudo_absences
cat("Тип объекта pseudo_absences:", class(pseudo_absences), "\n")
print(head(pseudo_absences))

# Если тип не SpatVector или есть проблемы, преобразуем заново
if (!inherits(pseudo_absences, "SpatVector")) {
  cat("Преобразование pseudo_absences в SpatVector...\n")
  pseudo_absences <- vect(as.data.frame(pseudo_absences), geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")
}

# Извлечение данных для псевдоотсутствий из всех слоев
env_pseudo_all <- terra::extract(env_layers_all, pseudo_absences)[, -1]  # Удаляем первый столбец (ID)

# Проверка на NA
if (any(is.na(env_pseudo_all))) {
  na_count <- sum(!complete.cases(env_pseudo_all))
  warning("Обнаружены NA в данных псевдоотсутствий. Удалено строк: ", na_count)
  env_pseudo_all <- na.omit(env_pseudo_all)
  pseudo_absences <- pseudo_absences[complete.cases(terra::extract(env_layers_all, pseudo_absences)[, -1]), ]
  pseudo_coords <- crds(pseudo_absences)
}
gc()

# Проверка результатов
cat("Успешно обработано", nrow(env_pseudo_all), "псевдоотсутствий\n")
print(head(env_pseudo_all))

# Сохранение псевдоотсутствий в CSV
pseudo_df <- data.frame(x = pseudo_coords[, 1], y = pseudo_coords[, 2])
write.csv(pseudo_df, "pseudo_absences.csv", row.names = FALSE)

# Сохранение псевдоотсутствий как shapefile
pseudo_sf <- st_as_sf(pseudo_absences)
st_write(pseudo_sf, "pseudo_absences.shp", delete_layer = TRUE)

# Визуализация
plot(env_layers_mask[[1]], main = "Study Area with Presence and Pseudo-Absences")
plot(polygon_vect, add = TRUE, border = "black")
plot(presence_vect, col = "red", pch = 16, cex = 0.8, add = TRUE)
if (exists("pseudo_absences")) {
  plot(pseudo_absences, col = "blue", pch = 16, cex = 0.8, add = TRUE)
} else {
  cat("Псевдоотсутствия не созданы из-за недостатка данных в маске.\n")
}

# Подготовка данных для анализа
# Создаем объединенный data.frame: присутствия (nf=1) и псевдоотсутствия (nf=0)
train_data <- as.data.frame(env_presence_all)
train_data$nf <- 1

pseudo_data <- as.data.frame(env_pseudo_all)
pseudo_data$nf <- 0
combined_data <- rbind(train_data, pseudo_data)
