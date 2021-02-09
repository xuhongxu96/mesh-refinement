#pragma once

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkInteractorStyleTerrain.h>
#include <vtkPolyData.h>

#include <unordered_map>

class MouseInteractorStyle : public vtkInteractorStyleTerrain {
 public:
  static MouseInteractorStyle* New();

  virtual void OnLeftButtonDown() override;

  std::unordered_map<vtkRenderer*, vtkPolyData*> DataMap;

  vtkNew<vtkDataSetMapper> selectedMapper;
  vtkNew<vtkActor> selectedActor;
};
