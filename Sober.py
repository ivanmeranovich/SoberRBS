#!/usr/bin/env python3
import argparse
import sys
import logging
import time
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def compute_similarity(seq1, seq2):
    """
    Вычисляет процент сходства двух последовательностей по посимвольному сравнению.
    Обе последовательности должны иметь одинаковую длину.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Последовательности должны иметь одинаковую длину!")
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100

def find_similar_occurrences(seq, motif, similarity_threshold=80.0):
    """
    Выполняет скользящий поиск по всей последовательности seq.
    Для каждого окна длины len(motif) вычисляется процент сходства с искомой последовательностью.
    Возвращает список кортежей (позиция, кандидат, сходство) для окон,
    где процент сходства не ниже similarity_threshold.
    """
    occurrences = []
    window_size = len(motif)
    for i in range(len(seq) - window_size + 1):
        candidate = seq[i:i+window_size]
        sim = compute_similarity(candidate, motif)
        if sim >= similarity_threshold:
            occurrences.append((i, candidate, sim))
    return occurrences

def is_valid_candidate(candidate_start, candidate_end, cds_intervals, threshold=40):
    """
    Проверяет, проходит ли кандидат следующие условия:
      1. Если кандидат полностью находится внутри любого CDS, он не подходит.
      2. Минимальное расстояние от кандидата до ближайшей границы CDS (начало или конец)
         должно быть не более threshold нуклеотидов.
         
    Возвращает кортеж (valid, min_distance),
      где valid – True, если кандидат допустим, а min_distance – минимальное расстояние до ближайшей границы CDS.
    """
    # Если в геноме нет CDS, условие считается невыполненным – кандидат не аннотируется.
    if not cds_intervals:
        return False, None

    # Если кандидат целиком находится внутри CDS, исключаем его.
    for (cds_start, cds_end) in cds_intervals:
        if candidate_start >= cds_start and candidate_end <= cds_end:
            return False, 0

    # Вычисляем минимальное расстояние до ближайшей границы CDS для каждого интервала.
    min_distance = None
    for (cds_start, cds_end) in cds_intervals:
        if candidate_end <= cds_start:
            dist = cds_start - candidate_end
        elif candidate_start >= cds_end:
            dist = candidate_start - cds_end
        else:
            # Этот случай не должен возникать, так как при пересечении кандидат уже исключён.
            dist = 0  
        if min_distance is None or dist < min_distance:
            min_distance = dist

    if min_distance is not None and min_distance <= threshold:
        return True, min_distance
    return False, min_distance

def main():
    parser = argparse.ArgumentParser(
        description="Скрипт ищет окно по заданной нуклеотидной последовательности (мотив) в GenBank-файле.\n"
                    "Аннотация происходит как фича типа RBS, если:\n"
                    "  • сходство с мотивом ≥ 80%,\n"
                    "  • (в обычном режиме) кандидат не находится внутри CDS,\n"
                    "     и расстояние до ближайшей границы CDS не превышает 40 нуклеотидов.\n\n"
                    "При использовании флага --test аннотируются все вхождения (с сходством ≥80%),\n"
                    "условия по CDS и дистанции игнорируются.\n\n"
                    "Флаг --hepl выводит расширенную справку по работе скрипта.\n\n"
                    "Пример вызова: python script.py input.gb AGGAGG output.gb [--test] [--hepl]"
    )
    parser.add_argument("input_file", help="Путь к входному GenBank-файлу")
    parser.add_argument("search_sequence", help="Искомая нуклеотидная последовательность (например, AGGAGG)")
    parser.add_argument("output_file", help="Путь к выходному GenBank-файлу с новыми аннотациями")
    parser.add_argument("--test", action="store_true",
                        help="Тестовый режим: аннотируются все вхождения (с сходством ≥80%), условия по CDS и дистанции игнорируются.")
    parser.add_argument("--hepl", action="store_true",
                        help="Вывод расширенной справки о работе скрипта.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Вывод подробной отладочной информации")
    args = parser.parse_args()

    if args.hepl:
        print("Расширенная справка:\n")
        print("Этот скрипт выполняет поиск окна с заданной нуклеотидной последовательностью (мотивом) по всему геному,\n"
              "вычисляет процент сходства найденного окна с искомой последовательностью и, если сходство не ниже 80%,\n"
              "проверяет, чтобы кандидат не находился внутри CDS и чтобы расстояние до ближайшей границы CDS не превышало 40 нуклеотидов.\n"
              "В обычном режиме аннотируются только кандидаты, удовлетворяющие этим условиям.\n"
              "При включении флага --test условия по CDS и расстоянию игнорируются и аннотируются все вхождения с нужным сходством.\n"
              "\nПример использования:\n"
              "  python script.py input.gb AGGAGG output.gb         (обычный режим)\n"
              "  python script.py input.gb AGGAGG output.gb --test    (тестовый режим)\n")
        sys.exit(0)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(message)s")
    start_time = time.time()

    try:
        record = SeqIO.read(args.input_file, "genbank")
    except Exception as e:
        sys.exit(f"Ошибка при чтении файла GenBank: {e}")

    genome_seq = str(record.seq).upper()
    motif = args.search_sequence.upper()
    similarity_threshold = 80.0
    proximity_threshold = 40  # порог изменен с 100 на 40 нуклеотидов

    logging.info(f"Входной файл: {args.input_file}")
    logging.info(f"Искомая последовательность: {motif}")
    logging.info(f"Порог сходства: {similarity_threshold}%")
    logging.info(f"Порог расстояния до CDS: {proximity_threshold} нт\n")

    # Получаем интервалы для всех CDS из аннотаций
    cds_intervals = []
    for feature in record.features:
        if feature.type.upper() == "CDS":
            cds_start = int(feature.location.start)
            cds_end = int(feature.location.end)
            cds_intervals.append((cds_start, cds_end))
    logging.info(f"Найдено CDS: {len(cds_intervals)}")

    # Поиск кандидатов по сходству
    candidates = find_similar_occurrences(genome_seq, motif, similarity_threshold=similarity_threshold)
    logging.info(f"Найдено кандидатов по сходству: {len(candidates)}\n")

    annotated_count = 0
    if args.test:
        logging.info("Режим --test: аннотируются все найденные вхождения (с сходством ≥80%), условия по CDS и дистанции игнорируются.")
        for pos, candidate_seq, sim in candidates:
            candidate_start = pos
            candidate_end = pos + len(motif)
            note = (f"TEST: Автоматически аннотированный RBS (иск. последовательность: {motif}, "
                    f"сходство: {sim:.2f}%)")
            new_feature = SeqFeature(
                FeatureLocation(candidate_start, candidate_end, strand=1),
                type="RBS",
                qualifiers={"note": note, "search_sequence": motif, "similarity": f"{sim:.2f}%"}
            )
            record.features.append(new_feature)
            annotated_count += 1
            logging.debug(f"TEST: Добавлена RBS: {candidate_start}-{candidate_end}; сходство: {sim:.2f}%")
    else:
        for pos, candidate_seq, sim in candidates:
            candidate_start = pos
            candidate_end = pos + len(motif)
            valid, distance = is_valid_candidate(candidate_start, candidate_end, cds_intervals, threshold=proximity_threshold)
            if not valid:
                logging.debug(f"Пропущено вхождение {candidate_start}-{candidate_end} (расстояние до CDS: {distance} нт)")
                continue
            note = (f"Автоматически аннотированный RBS (иск. последовательность: {motif}, сходство: {sim:.2f}%)"
                    f"; расстояние до ближайшей границы CDS: {distance} нт")
            new_feature = SeqFeature(
                FeatureLocation(candidate_start, candidate_end, strand=1),
                type="RBS",
                qualifiers={"note": note, "search_sequence": motif, "similarity": f"{sim:.2f}%", "distance": str(distance)}
            )
            record.features.append(new_feature)
            annotated_count += 1
            logging.debug(f"Добавлена RBS: {candidate_start}-{candidate_end}; сходство: {sim:.2f}%; расстояние: {distance} нт")

    logging.info(f"\nАннотировано вхождений RBS: {annotated_count}")

    try:
        SeqIO.write(record, args.output_file, "genbank")
        logging.info(f"\nРезультаты записаны в файл: {args.output_file}")
    except Exception as e:
        sys.exit(f"Ошибка при записи файла GenBank: {e}")

    end_time = time.time()
    logging.info(f"Время выполнения: {end_time - start_time:.2f} секунд")

if __name__ == "__main__":
    main()
