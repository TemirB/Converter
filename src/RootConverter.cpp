#include "RootConverter.h"
#include "McDstWriter.h"

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <vector>

bool ConvertRootToMcDst(const std::string& inName, const std::string& outName) {
    try {
        if (std::filesystem::weakly_canonical(inName) ==
            std::filesystem::weakly_canonical(outName)) {
            throw std::runtime_error("input and output must be different files");
        }
        std::unique_ptr<TFile> input(TFile::Open(inName.c_str(), "READ"));
        if (!input || input->IsZombie())
            throw std::runtime_error("cannot open input ROOT file");
        auto* tree = dynamic_cast<TTree*>(input->Get("Particles"));
        if (!tree)
            throw std::runtime_error("missing Particles tree");
        for (const auto* name : {"npart", "impactParameter", "t", "x", "y", "z", "mass", "p0", "px",
                                 "py", "pz", "pdg", "ID", "charge"}) {
            if (!tree->GetBranch(name))
                throw std::runtime_error("missing ROOT branch: " + std::string(name));
        }
        TTreeReader reader(tree);
        TTreeReaderValue<int> count(reader, "npart");
        TTreeReaderValue<double> impact(reader, "impactParameter");
        TTreeReaderValue<std::vector<double>> t(reader, "t"), x(reader, "x"), y(reader, "y"),
            z(reader, "z");
        TTreeReaderValue<std::vector<double>> mass(reader, "mass"), energy(reader, "p0");
        TTreeReaderValue<std::vector<double>> px(reader, "px"), py(reader, "py"), pz(reader, "pz");
        TTreeReaderValue<std::vector<int>> pdg(reader, "pdg"), id(reader, "ID"),
            charge(reader, "charge");
        McDstWriter writer(outName, "ROOT final particle lists; Npart/Ncoll and spectator status "
                                    "unknown; coordinates from input; genealogy unknown");
        if (!writer.Error().empty())
            throw std::runtime_error(writer.Error());
        int number = 0;
        while (reader.Next()) {
            if (count.GetSetupStatus() < 0 || impact.GetSetupStatus() < 0 ||
                t.GetSetupStatus() < 0 || x.GetSetupStatus() < 0 || y.GetSetupStatus() < 0 ||
                z.GetSetupStatus() < 0 || mass.GetSetupStatus() < 0 ||
                energy.GetSetupStatus() < 0 || px.GetSetupStatus() < 0 || py.GetSetupStatus() < 0 ||
                pz.GetSetupStatus() < 0 || pdg.GetSetupStatus() < 0 || id.GetSetupStatus() < 0 ||
                charge.GetSetupStatus() < 0)
                throw std::runtime_error("invalid ROOT input branches");
            if (*count < 0 || !std::isfinite(*impact) || *impact < 0.)
                throw std::runtime_error("invalid ROOT event metadata");
            const auto n = static_cast<std::size_t>(*count);
            for (const auto size :
                 {t->size(), x->size(), y->size(), z->size(), mass->size(), energy->size(),
                  px->size(), py->size(), pz->size(), pdg->size(), id->size(), charge->size()}) {
                if (size != n)
                    throw std::runtime_error("particle vector size differs from npart");
            }
            OscarEvent event;
            event.number = number++;
            event.impact_parameter = *impact;
            event.initial_nucleons = event.participants = event.spectators = -1;
            event.binary_collisions = event.nn_elastic_collisions = event.nn_inelastic_collisions =
                -1;
            std::unordered_set<int> seen_ids;
            for (std::size_t i = 0; i < n; ++i) {
                Particle p;
                p.t = (*t)[i];
                p.x = (*x)[i];
                p.y = (*y)[i];
                p.z = (*z)[i];
                p.mass = (*mass)[i];
                p.p0 = (*energy)[i];
                p.px = (*px)[i];
                p.py = (*py)[i];
                p.pz = (*pz)[i];
                p.pdg = (*pdg)[i];
                p.id = (*id)[i];
                p.charge = (*charge)[i];
                for (double value : {p.t, p.x, p.y, p.z, p.mass, p.p0, p.px, p.py, p.pz}) {
                    if (!std::isfinite(value))
                        throw std::runtime_error("non-finite ROOT particle value");
                }
                if (p.id < 0 || p.pdg == 0 || p.mass < 0. || p.p0 <= 0.)
                    throw std::runtime_error("invalid ROOT particle record");
                if (!seen_ids.insert(p.id).second)
                    throw std::runtime_error("duplicate ROOT particle ID");
                event.final_time = std::max(event.final_time, p.t);
                event.particles.push_back(p);
            }
            if (!writer.WriteEvent(event))
                throw std::runtime_error(writer.Error());
        }
        if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && tree->GetEntries() != 0)
            throw std::runtime_error("ROOT input read failed");
        if (!writer.Finish(number, 0))
            throw std::runtime_error(writer.Error());
        std::cout << "Converted " << number << " ROOT events.\n";
        return true;
    } catch (const std::exception& error) {
        std::cerr << "Conversion failed: " << error.what() << '\n';
        return false;
    }
}
