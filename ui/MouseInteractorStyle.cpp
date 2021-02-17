#include "MouseInteractorStyle.h"

#include <vtkCellPicker.h>
#include <vtkExtractSelection.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkUnstructuredGrid.h>

void MouseInteractorStyle::OnLeftButtonDown() {
  // Get the location of the click (in window coordinates)
  int* pos = this->GetInteractor()->GetEventPosition();
  this->FindPokedRenderer(pos[0], pos[1]);

  vtkNew<vtkCellPicker> picker;
  picker->SetTolerance(0.0005);

  // Pick from this location.
  picker->Pick(pos[0], pos[1], 0, this->GetCurrentRenderer());

  double* worldPosition = picker->GetPickPosition();
  std::cout << "Cell id is: " << picker->GetCellId() << std::endl;

  auto mesh = this->DataMap[this->GetCurrentRenderer()];

  if (picker->GetCellId() != -1) {
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
    extractSelection->SetInputData(0, mesh);
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());

    std::cout << "There are " << selected->GetNumberOfPoints()
              << " points in the selection." << std::endl;

    std::cout.precision(8);
    std::cout.setf(std::ios::fixed);

    // Output angles
    double p[3][3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        std::cout << (p[i][j] = selected->GetPoint(i)[j]) << " ";
    std::cout << std::endl;

    double vec[3][3];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 2; ++j) {
        vec[i][j] = p[i][j] - p[(i + 1) % 3][j];
      }
      vec[i][2] = 0;  // calculate angle in 2d only
    }

    for (int i = 0; i < 3; ++i) {
      double current_radian = vtkMath::Pi() - vtkMath::AngleBetweenVectors(
                                                  vec[i], vec[(i + 1) % 3]);
      std::cout << vtkMath::DegreesFromRadians(current_radian) << std::endl;
    }

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

vtkStandardNewMacro(MouseInteractorStyle);
