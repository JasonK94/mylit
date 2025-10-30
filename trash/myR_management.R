library(devtools)

# 패키지 로드
load_all("myR")

# myR 패키지만 빌드
build(pkg = "myR")   

# myR 설치
install(pkg = "myR")

# myR R CMD check
check(pkg = "myR")

# myR 패키지 제거
remove.packages("myR")

# myR 패키지 설치 확인
installed.packages()
