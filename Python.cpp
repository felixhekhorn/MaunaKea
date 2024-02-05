#include <pybind11/pybind11.h>

#include "./MaunaKea.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MaunaKea, m) {
  m.doc() = "Hadroproduction of heavy quark flavors";

  py::class_<MaunaKea::IntegrationConfig>(m, "IntegrationConfig")
      .def_readwrite("verbosity", &MaunaKea::IntegrationConfig::verbosity, py::doc("level of output"))
      .def_readwrite("calls", &MaunaKea::IntegrationConfig::calls, py::doc("calls"))
      .def_readwrite("warmupCalls", &MaunaKea::IntegrationConfig::warmupCalls, py::doc("calls for warmup"))
      .def_readwrite("iterations", &MaunaKea::IntegrationConfig::iterations, py::doc("iterations"))
      .def_readwrite("warmupIterations", &MaunaKea::IntegrationConfig::warmupIterations,
                     py::doc("iterations during warmup"))
      .def_readwrite("adaptChi2", &MaunaKea::IntegrationConfig::adaptChi2, py::doc("iterate until |chi2-1| < 0.5?"))
      .def_readwrite("bins", &MaunaKea::IntegrationConfig::bins, py::doc("number of bins"))
      .def("toString", &MaunaKea::IntegrationConfig::toString, py::doc("Dump current state to string"))
      .doc() = "Integration parameter configuration";

  py::class_<MaunaKea::IntegrationOutput>(m, "IntegrationOutput")
      .def_readonly("result", &MaunaKea::IntegrationOutput::result, py::doc("result"))
      .def_readonly("error", &MaunaKea::IntegrationOutput::error, py::doc("(absolute) error"))
      .def_readonly("chi2", &MaunaKea::IntegrationOutput::chi2, py::doc("chi^2"))
      .def_readonly("chi2iter", &MaunaKea::IntegrationOutput::chi2iter,
                    py::doc("number of iteration to converge chi^2"))
      .def("toString", &MaunaKea::IntegrationOutput::toString, py::doc("Dump current state to string"))
      .doc() = "Integration result";

  py::class_<MaunaKea::MaunaKea>(m, "MaunaKea")
      .def(py::init<cdbl, cuint, cuint, cuint>())
      .def_readwrite("intCfg", &MaunaKea::MaunaKea::intCfg, py::doc("Integration configuration"))
      .def("setHadronicS", &MaunaKea::MaunaKea::setHadronicS, py::doc("Set hadronic Mandelstam S_h"))
      .def("setScaleRatios", &MaunaKea::MaunaKea::setScaleRatios,
           py::doc("Set renormalization and factorization scale ratios xi_{R/F} = mu_{R/F}/m"))
      .def("setGridCentralScaleRatio", &MaunaKea::MaunaKea::setGridCentralScaleRatio,
           py::doc("Set grid central scale ratio xi = mu/m"))
      .def("setPDF", &MaunaKea::MaunaKea::setPDF, py::doc("Set reference PDF"))
      .def("run", &MaunaKea::MaunaKea::run, py::doc("Run calculation"))
      .def("write", &MaunaKea::MaunaKea::write, py::doc("Write grid to disk"))
      .def("getIntegrationOutput", &MaunaKea::MaunaKea::getIntegrationOutput,
           py::doc("Copy of current integration output"), py::return_value_policy::copy)
      .doc() = "Main application class";
}
