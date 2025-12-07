#include<bits/stdc++.h>
using namespace std;

class AbstarctEmployee {
    virtual void AskForPromotion() = 0;
};

class Employee : AbstarctEmployee {
private:
    string Name;
    string Company;
    int Age;
public:
    void setName(string name){
        Name = name;
    }
    void setCompany(string company){
        Company = company;
    }
    void setAge(int age){
        Age = age;
    }
    string getName(){
        return Name;
    }
    string getCompany(){
        return Company;
    }
    int getAge(){
        return Age;
    }
    void IntroduceYoureself(){
        cout << "Name-" << Name << "\n";
        cout << "Company-"<< Company << "\n";
        cout << "Age-"<< Age << "\n";
    }
    Employee(string name, string company, int age){
        Name = name;
        Company = company;
        Age = age;
    }
    void AskForPromotion();
};

void Employee::AskForPromotion() {
        if (Age > 24){
            cout << Name << " got promoted!" << "\n";
        }
        else{
            cout << Name << ", sorry NO promotion!" << "\n";
        }
    }

class Developer: Employee {
public:
    string FavProgrammingLanguage;
    Developer(string name, string company, int age, string favProgrammingLanguage)
        :Employee(name, company, age)
    {
        FavProgrammingLanguage = favProgrammingLanguage;
    }
    void FixBug() {
        cout << getName() << " fix bug using " << FavProgrammingLanguage << "\n";
    }
};

int main(){
    Developer* d = new Developer("Ankit", "xyz", 23, "c++");
    d->FixBug();
    Employee* employee1 = new Employee("Ankit", "xyz", 23);
    Employee* employee2 = new Employee("Raunak", "yzx", 25);
    employee1->AskForPromotion();
    employee2->AskForPromotion();
}