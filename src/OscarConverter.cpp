#include "OscarConverter.h"
#include "Particle.h"
#include "Mode.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "McRun.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McArrays.h"

// Кол-во начальных частиц
// Далее по нему будет определяться является ли частица начальной
int kInitNucl = 394;

// Хелпер для записи события
void processEvent(
    TClonesArray* arrays[McArrays::NAllMcArrays],
    const std::unordered_map<int, Particle>& out,
    double impactParameter, int evNum, int nPart
) {
    // [Future] Можно из interaction поддерживать сбор child
    int child[2] = {-1, 1};
    int idx = 0;

    // пред запрос
    auto* particleArray = arrays[McArrays::Particle];

    for (const auto & entry : out) {
        int entryIdx = particleArray->GetEntries();
        const Particle& part = entry.second;
        // Use placement new: TClonesArray pre-allocates memory,
        // and we construct McParticle directly in the next free slot.
        // arrays[...]->GetEntries() gives the index of that slot,
        // and operator[] returns a void* to the raw memory location.
        
        new((*particleArray)[entryIdx]) McParticle(
            idx++, part.pdg, 0, 0, 0, -1, 0, child,
            part.px, part.py, part.pz, part.p0,
            part.x, part.y, part.z, part.t
        );
    }

    McEvent* event = new((*(arrays[McArrays::Event]))[arrays[McArrays::Event]->GetEntries()]) McEvent();

    event->setEventNr(evNum);
    event->setB(impactParameter);
    event->setNpart(nPart);
    event->setPhi(0.);
    event->setNes(1);
    event->setComment("");
    event->setStepNr(1);
    event->setStepT(200.);
}

// Перписаный конвертер:
// С прошлой версии поменял условие отбора спектатора.
// А также перенес обработку out-end блока в разбор end строчки,
// чтобы логически отделить код записи в буферы и записи в 
bool OscarConverter::Convert(const std::string& nInput, const std::string& nOutput) {
    // Создание входного файла и проверка, что он создался
    TFile* output = TFile::Open(nOutput.c_str(), "RECREATE", "Oscar to McDst");
    if (!output || output->IsZombie()) {
        std::cerr << "Cannot create output file: " << nOutput << std::endl;
        return false;
    }
    
    // Вроде как эта строка оптимизирует сжимаемость
    // Должно быть полезно, при обработке большого кол-ва файлов
    output->SetCompressionLevel(1);

    // Создаем McDst дерево
    TTree* tree = new TTree("McDst", "Converted Oscar Data");
    TClonesArray* arrays[McArrays::NAllMcArrays];

    // Иннициализируем базовые массивы
    for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) {
        arrays[i] = new TClonesArray(McArrays::mcArrayTypes[i], McArrays::mcArraySizes[i]);
        arrays[i]->SetOwner(kFALSE);
        auto br = tree->Branch(McArrays::mcArrayNames[i], &arrays[i], 65536, 99);
        if (br) br->SetAutoDelete(kFALSE);
    }

    // Открываем входной файл
    std::ifstream input(nInput);
    if (!input.is_open()) {
        std::cerr << "Cannot open input file: " << nInput << std::endl;
        return false;
    }

    // Общие параметры события
    // - isElastic = флаг события, определяющий является ли оно упругим
    // - evNum = номер события
    // - nPart = кол-во частиц в событии, используется для проверки события на упругость
    bool isElastic = false;

    int evNum = -1;
    int nPart = -1;

    // [Future] Параметры взаимодействия
    // - iIn = кол-во входных частиц
    // - iOut = кол-во выходных частиц 
    // int iIn = -1;
    // int iOut = -1;

    // Буферы
    // - interaction = буфер частиц из взаимодействий
    //                 отслеживает все уникальные частицы из каждого взаимодействия
    //                 отчистка происходит в момент завершения всего события
    // - out = буфер частиц, доживших до 200 fm/c
    //         основной буфер для записи в дерево
    //         отчистка происходит в момент завершения всего события
    std::unordered_map<int, Particle> interactions;
    std::unordered_map<int, Particle> out;
    // Резервируем память, оценочное значение.
    interactions.reserve(15000);
    out.reserve(15000);

    // Основной принцип работы программы
    // - mode = содержит в себе текущий режим работы программы,
    //          меняется в зависимости от ключевого слова:
    //          - # interaction -> Mode::Interaction 
    //              Режим записи в буфер interation
    //          - # ... in -> Mode::InEvent
    //              Вероятно лишний режим работы, но влияния на корректность и производительность не должно быть, 
    //              пока оставляю, после валидации нужно будет убрать
    //          - # ... out -> Mode::OutEvent 
    //              Режим записи частиц в буфер out, поддерживая уникальность каждого ID, 
    //              посредством записи только последней частицы с этим ID,
    //              а также проверяем, является ли событие упругим
    //          - # ... end -> Mode::EndEvent 
    //              По сути тоже лишний режим (перенос обработки в парсинг облегчает логику)
    //              Режим обработки, проставляем в out спектаторы и записываем частицы из out в McDst массивы.
    Mode mode = Mode::Init;
    
    std::string line;
    line.reserve(256);

    while (std::getline(input, line)) {
        if (line.empty()) continue;

        // std::istringstream iss(line);

        // Разбираем по словам комментарий:
        if (line[0] == '#') {
            char word1[32], word2[32];
            // std::string interaction, dummy, keyWord;

            if (sscanf(line.c_str(), "# %31s %31s", word1, word2) < 2) {
                continue; // Не хватает слов, пропускаем
            }
            //       #      event
            //       #   interaction
            // if (dummy == "interaction") {
            if (strcmp(word1, "interaction") == 0) {
                // Поиск кол-ва частиц в взаимодействии
                // [Future] можно будет поддерживать child, etc
                // iss >> dummy >> iIn >> dummy >> iOut;
                //      in       2      out      2  ... (rho    0.0000000 weight     39.71011 partial   39.7101111 type     1)
                mode = Mode::Interaction;
            // } else if (dummy == "event") {
            } else if (strcmp(word1, "event") == 0) {
                char keyword[4];
                // Поиск ключевых слов
                if (sscanf(line.c_str(), "# event %d %3s %d", &evNum, keyword, &nPart) != 3) {
                    mode = Mode::SkipEvent;
                    continue;
                }
                // iss >> evNum >> keyWord >> nPart;
                //       0        in        394
                //       0        out       421
                //       0        end        0
                // if (keyWord == "in") {
                if (strcmp(keyword, "in") == 0) {
                    mode = Mode::InEvent;
                    // [OUTDATED] Не нужный режим работы
                // } else if (keyWord == "out") {
                } else if (strcmp(keyword, "out") == 0) {
                    mode = (nPart == kInitNucl) ? Mode::SkipEvent : Mode::OutEvent;
                    isElastic = (mode == Mode::SkipEvent);
                // } else if (keyWord == "end") {
                } else if (strcmp(keyword, "end") == 0) {
                    mode = Mode::EndEvent;
                    double impactParameter = -1.;
                    char yesno[4];

                    if (sscanf(line.c_str(), "# event %*d end %*d impact %lf scattering_projectile_target %3s", &impactParameter, yesno) == 2) {
                        if (isElastic && strcmp(yesno, "no") == 0) {
                            impactParameter = -1.;
                        }
                    }
                    // iss >> dummy >> impactParameter >>             dummy       >> dummy; 
                    //    impact        11.095         scattering_projectile_target yes

                    // scattering_projectile_target - показывает упругость события,
                    // yes - событие было, новые частицы родились
                    // no - небыло новых частиц кроме начальных = упругое событие
                    // if (isElastic && dummy == "no") {
                    //     impactParameter = -1.;
                    // }

                    if (!isElastic) {
                        processEvent(arrays, out, impactParameter, evNum, nPart);
                        tree->Fill();
                    }

                    // clear
                    interactions.clear();
                    out.clear();
                    for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) arrays[i]->Clear();
                    isElastic = false;
                } else {
                    mode = Mode::SkipEvent;
                }
            }
            // Распарсили строку, переходим дальше
            // Или это сервисный комментарий
            continue;
        }

        if (mode == Mode::SkipEvent) continue;

        Particle p;
        if (sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
                   &p.t, &p.x, &p.y, &p.z, 
                   &p.mass, &p.p0, &p.px, &p.py, &p.pz,
                   &p.pdg, &p.ID, &p.charge) != 12) {
            continue;
        }

        switch (mode) {
            case Mode::Interaction:
                interactions.try_emplace(p.ID, p);
                break;

            case Mode::OutEvent: {
                if (isElastic) continue;

                auto it = interactions.find(p.ID);
                bool has_interaction = (it != interactions.end());
                
                if (p.ID < kInitNucl && !has_interaction) {
                    p.SetSpectator();
                } else if (has_interaction) {
                    p = it->second;
                }
                out.try_emplace(p.ID, std::move(p));
                
                break;
            }

            default:
                break;
        }
    }


    McRun run("SMASH", "Converted from Oscar file",
              0, 0, 0., 0, 0, 0.,
              0., 0., -1, 0, 0, 0., tree->GetEntries());
    run.Write();

    std::cout << "Conversion completed. Total events: " << tree->GetEntries() << std::endl;

    output->Write();
    output->Close();

    return true;
}


// // [OUTDATED]
// bool oldConvert(const std::string& inputFilename, const std::string& outputFilename) {
//     TFile* outputFile = TFile::Open(outputFilename.c_str(), "RECREATE", "Oscar to McDst");
//     if (!outputFile || outputFile->IsZombie()) {
//         std::cerr << "Cannot create output file: " << outputFilename << std::endl;
//         return false;
//     }
//     outputFile->SetCompressionLevel(1);

//     TTree* tree = new TTree("McDst", "Converted Oscar Data");
//     TClonesArray* arrays[McArrays::NAllMcArrays];

//     for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) {
//         arrays[i] = new TClonesArray(McArrays::mcArrayTypes[i], McArrays::mcArraySizes[i]);
//         arrays[i]->SetOwner(kFALSE);
//         auto br = tree->Branch(McArrays::mcArrayNames[i], &arrays[i], 65536, 99);
//         if (br) br->SetAutoDelete(kFALSE);
//     }

//     std::ifstream infile(inputFilename);
//     if (!infile.is_open()) {
//         std::cerr << "Cannot open input file: " << inputFilename << std::endl;
//         return false;
//     }

//     bool isElastic = false;
//     int ev_num = -1;
//     double n_part = -1;
//     Mode mode = Mode::Init;
//     int startParticlesNum = 394;

//     double timpactParameter = -1.;

//     std::unordered_map<int, Particle> buffer;
//     std::unordered_map<int, Particle> endBuffer;
//     std::unordered_map<int, Particle> eventBuffer;

//     std::string line;
//     while (std::getline(infile, line)) {
//         if (line.empty()) continue;

//         std::istringstream iss(line);

//         if (line[0] == '#') {
//             std::string interaction, dummy, keyWord;
//             iss >> dummy >> dummy >> ev_num >> keyWord >> n_part;
//             //       #      event       0       in          394
//             if (dummy == "interaction") {
//                 mode = Mode::Interaction;
//                 continue;
//             } else if (dummy == "event") {
//                 if (keyWord == "in") {
//                     mode = Mode::InEvent;
//                     isElastic = false;

//                     buffer.clear();
//                     endBuffer.clear();
//                     eventBuffer.clear();
//                     for (unsigned int i = 0; i < McArrays::NAllMcArrays; i++) arrays[i]->Clear();
//                     continue;
//                 } else if (keyWord == "out") {
//                     mode = (int(n_part) == startParticlesNum) ? Mode::SkipEvent : Mode::OutEvent;
//                     isElastic = (mode == Mode::SkipEvent);
//                     continue;
//                 } else if (keyWord == "end") {
//                     mode = Mode::EndEvent;
//                     timpactParameter = n_part;
//                     if (isElastic) {
//                         timpactParameter = -1.;
//                         continue;
//                     }
//                 } else {
//                     mode = Mode::SkipEvent;
//                     continue;
//                 }
//             }
//         }

//         if (mode == Mode::SkipEvent) continue;

//         double t, x, y, z, mass, p0, px, py, pz;
//         int pdg, ID, charge;
//         iss >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pdg >> ID >> charge;
//         Particle p(t, x, y, z, mass, p0, px, py, pz, pdg, ID, charge);

//         if (mode == Mode::Interaction || mode == Mode::InEvent) {
//             if (ID < 394) p.participant = true;
//             buffer[ID] = p;
//         } else if (mode == Mode::OutEvent) {
//             if (isElastic) continue;
//             endBuffer[ID] = p;
//         } else if (mode == Mode::EndEvent) {
//             for (const auto& entry : buffer) {
//                 int id = entry.first;
//                 const Particle& initParticle = entry.second;

//                 if (endBuffer.find(id) == endBuffer.end()) continue;

//                 if (!initParticle.participant && initParticle.ID < 394) {
//                     Particle specParticle = initParticle;
//                     specParticle.SetSpectator();
//                     eventBuffer[id] = specParticle;
//                 } else {
//                     eventBuffer[id] = initParticle;
//                 }
//             }

//             int child[2] = {-1, -1};
//             int idx = 0;
//             for (const auto& entry : eventBuffer) {
//                 const Particle& p1 = entry.second;
//                 new((*(arrays[McArrays::Particle]))[arrays[McArrays::Particle]->GetEntries()])
//                 McParticle(idx++, p1.pdg, 0, 0, 0, -1, 0, child,
//                           p1.px, p1.py, p1.pz, p1.p0,
//                           p1.x, p1.y, p1.z, p1.t);
//             }

//             McEvent* event = new((*(arrays[McArrays::Event]))[arrays[McArrays::Event]->GetEntries()]) McEvent();
//             event->setEventNr(ev_num);
//             event->setB(timpactParameter);
//             event->setPhi(0.);
//             event->setNes(1);
//             event->setComment("");
//             event->setStepNr(1);
//             event->setStepT(200.);

//             tree->Fill();
//             buffer.clear();
//             endBuffer.clear();
//             eventBuffer.clear();
//         }
//     }

//     McRun run("SMASH", "Converted from Oscar file",
//               0, 0, 0., 0, 0, 0.,
//               0., 0., -1, 0, 0, 0., tree->GetEntries());
//     run.Write();

//     std::cout << "Conversion completed. Total events: " << tree->GetEntries() << std::endl;

//     outputFile->Write();
//     outputFile->Close();

//     return true;
// }
