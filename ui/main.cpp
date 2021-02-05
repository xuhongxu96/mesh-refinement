#include <MeshReader.h>
#include <Refiner.h>
#include <vtkActor.h>
#include <vtkCellPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTerrain.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSimplePointsReader.h>
#include <vtkUnstructuredGrid.h>

class MouseInteractorStyle : public vtkInteractorStyleTerrain {
 public:
  static MouseInteractorStyle* New();

  virtual void OnLeftButtonDown() override {
    // Get the location of the click (in window coordinates)
    int* pos = this->GetInteractor()->GetEventPosition();

    vtkNew<vtkCellPicker> picker;
    picker->SetTolerance(0.0005);

    // Pick from this location.
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    double* worldPosition = picker->GetPickPosition();
    std::cout << "Cell id is: " << picker->GetCellId() << std::endl;

    if (picker->GetCellId() != -1) {
      std::cout << "Pick position is: " << worldPosition[0] << " "
                << worldPosition[1] << " " << worldPosition[2] << endl;

      vtkNew<vtkIdTypeArray> ids;

      ids->SetNumberOfComponents(1);
      ids->InsertNextValue(picker->GetCellId());

      vtkSmartPointer<vtkSelectionNode> selectionNode =
          vtkSmartPointer<vtkSelectionNode>::New();
      selectionNode->SetFieldType(vtkSelectionNode::CELL);
      selectionNode->SetContentType(vtkSelectionNode::INDICES);
      selectionNode->SetSelectionList(ids);

      vtkSmartPointer<vtkSelection> selection =
          vtkSmartPointer<vtkSelection>::New();
      selection->AddNode(selectionNode);

      vtkSmartPointer<vtkExtractSelection> extractSelection =
          vtkSmartPointer<vtkExtractSelection>::New();
      extractSelection->SetInputData(0, this->Data);
      extractSelection->SetInputData(1, selection);
      extractSelection->Update();

      // In selection
      vtkSmartPointer<vtkUnstructuredGrid> selected =
          vtkSmartPointer<vtkUnstructuredGrid>::New();
      selected->ShallowCopy(extractSelection->GetOutput());

      std::cout << "There are " << selected->GetNumberOfPoints()
                << " points in the selection." << std::endl;
      std::cout << "There are " << selected->GetNumberOfCells()
                << " cells in the selection." << std::endl;
      selectedMapper->SetInputData(selected);
      selectedActor->SetMapper(selectedMapper);
      selectedActor->GetProperty()->EdgeVisibilityOn();
      selectedActor->GetProperty()->SetColor(0, 200, 0);

      selectedActor->GetProperty()->SetLineWidth(3);

      this->Interactor->GetRenderWindow()
          ->GetRenderers()
          ->GetFirstRenderer()
          ->AddActor(selectedActor);
    }

    // Forward events
    vtkInteractorStyleTerrain::OnLeftButtonDown();
  }

  vtkPolyData* Data;
  vtkNew<vtkDataSetMapper> selectedMapper;
  vtkNew<vtkActor> selectedActor;
};

vtkStandardNewMacro(MouseInteractorStyle);

int main() {
  vtkNew<vtkSimplePointsReader> point_reader;
  point_reader->SetFileName(DATA_PATH "vectors.txt");
  point_reader->Update();

  mr::MeshReader mesh_reader;
  mr::Mesh mesh = mesh_reader.Read(DATA_PATH "original.mesh");

  mr::Refiner refiner(point_reader->GetOutput());
  auto own_colors = refiner.Refine(mesh.poly_data);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(mesh.poly_data);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetSpecular(.4);
  actor->GetProperty()->SetSpecularPower(30.0);

  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Mesh Refinement");

  vtkNew<vtkInteractorStyleTerrain> style;
  // style->SetDefaultRenderer(renderer);
  // style->Data = mesh.poly_data;

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetInteractorStyle(style);

  renderer->AddActor(actor);
  renderer->SetBackground(255, 255, 255);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}