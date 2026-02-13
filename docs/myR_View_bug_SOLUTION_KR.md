# myR View() 버그 - 최종 심층 분석 및 해결

## 🔍 핵심 원인: `import(stats)`의 나비 효과

사용자님의 관찰과 테스트 결과, 문제는 **`NAMESPACE`에 `import(stats)`가 포함되었을 때** 발생했습니다.

### 1. 왜 `load_all()`에서만 문제가 발생했나?

`devtools::load_all()`은 설치된 패키지를 로드하는 `library()`와 다르게 동작합니다:
- **`library()`**: 컴파일되고 정리된 `NAMESPACE`를 정해진 순서대로 깔끔하게 부착합니다.
- **`load_all()`**: 소스 코드를 `source()`하여 개발용 임시 환경(shim)에 올립니다. 이 과정에서 의존성 패키지들(`Imports`)을 로드하는데, 순서와 우선순위가 미묘하게 달라질 수 있습니다.

### 2. `import(stats)`가 범인인 이유

실패했던 시점(`0e844ff`)의 `NAMESPACE`에는 `import(stats)`가 있었습니다.
이것이 다음과 같은 연쇄 반응을 일으켰습니다:

1. `stats` 패키지가 명시적으로 import 되면서, **`stats::filter`, `stats::lag`가 `dplyr`의 함수들을 덮어썼습니다(Overrided)**.
2. 이로 인해 패키지 로딩 및 함수 탐색 경로(Search Path)에 변동이 생겼습니다.
3. 이 틈을 타 `BiocGenerics` 패키지(Seurat 등의 의존성)가 정의한 **S4 Generic `as.data.frame`**이 우선권을 잡게 되었습니다.
4. 결과적으로 `View(seurat_obj)` 호출 시, `base::as.data.frame` 대신 `BiocGenerics`의 S4 제네릭이 호출되었고, Seurat 객체는 이를 처리하지 못해 에러가 발생했습니다.

### 3. `conflicted`가 막지 못한 이유

`conflicted` 패키지는 사용자가 함수를 호출할 때(`filter()`) 발생하는 충돌을 해결해줍니다. 하지만 이번 문제는 **패키지 로딩 시점에 시스템 내부적으로 `as.data.frame` 메서드 테이블이 꼬이는 문제**여서, `conflicted`가 개입하기도 전에 이미 상황이 종료되었을 가능성이 큽니다.

## ✅ 최종 해결책

**`import(stats)`를 제거하고, `import(dplyr)`는 유지합니다.**

- 현재 코드 베이스에는 `@import stats` 구문이 없으므로, `roxygen2::roxygenise()`를 실행하면 자동으로 `import(stats)`가 없는 깨끗한 `NAMESPACE`가 생성됩니다.
- 이것은 `View()`가 정상 작동했던 과거 시점(`4754f14`)과 동일한 구성입니다.
- `zzz.R`의 `conflicted` 설정은 안전장치(보험)로 유지합니다.

이제 `dplyr`의 모든 기능을 문제없이 사용하면서 `View()`도 정상 작동할 것입니다.
