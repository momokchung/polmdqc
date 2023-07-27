/////////////////////////////////////////////////////////
//                                                     //
//  getnumb.cpp  --  extract an integer from a string  //
//                                                     //
/////////////////////////////////////////////////////////

// "getnumb" searches an input string from left to right for an
// integer and puts the numeric value in "number"; returns zero
// with "next" unchanged if no integer value is found
//
// variables and parameters:
//
// string    input character string to be searched
// number    output with the first integer in the string
// next      input with first position of search string;
//             output with the position following the number


#include "getnumb.h"
#include <sstream>

void getnumb(std::string& string, int& number, int& next)
{
    number = 0;
    std::string newString = string.substr(next);
    std::istringstream iss(newString);

    if (!(iss >> number)) return;

    std::streampos pos = iss.tellg();
    int posInt = static_cast<int>(pos);
    int length = newString.length();
    int addLength = (posInt == -1 ? length : posInt);

    next += addLength;
}


//       subroutine getnumb (string,number,next)
//       use ascii
//       implicit none
//       integer i,j,length
//       integer number,digit
//       integer next,trimtext
//       integer first,last,code
//       integer initial,final
//       integer place(10)
//       logical positive
//       logical negative
//       logical numeral
//       character*1 letter
//       character*(*) string
//       data place  / 1, 10, 100, 1000, 10000, 100000, 1000000,
//      &              10000000, 100000000, 1000000000 /
// c
// c
// c     initialize number and get the input text string length
// c
//       number = 0
//       positive = .false.
//       negative = .false.
//       numeral = .false.
//       length = trimtext(string(next:))
// c
// c     search the string for the first run of numeric characters
// c
//       first = next
//       last = 0
//       initial = next
//       final = next + length - 1
//       do i = initial, final
//          letter = string(i:i)
//          code = ichar(letter)
//          if (letter.ge.'0' .and. letter.le.'9') then
//             if (.not. numeral) then
//                numeral = .true.
//                first = i
//             end if
//             if (i .eq. final) then
//                last = final
//                next = i + 1
//             end if
//          else if (code.eq.plus .and. .not.positive) then
//             positive = .true.
//          else if (code.eq.minus .and. .not.negative) then
//             negative = .true.
//          else if (numeral) then
//             if (code.eq.space .or. code.eq.tab .or.
//      &          code.eq.comma .or. code.eq.semicolon .or.
//      &          code.eq.colon .or. code.eq.underbar) then
//                last = i - 1
//                next = i
//             else
//                numeral = .false.
//             end if
//             goto 10
//          else if (positive .or. negative) then
//             numeral = .false.
//             goto 10
//          else if (code.ne.space .and. code.ne.tab) then
//             numeral = .false.
//             goto 10
//          end if
//       end do
//    10 continue
// c
// c     trim the actual number if it is too big to return
// c
//       if (.not. numeral)  next = initial
//       last = min(last,first+9)
// c
// c     convert the text numeral into an integer number
// c
//       j = 0
//       do i = last, first, -1
//          j = j + 1
//          if (string(i:i) .eq. '0') then
//             digit = 0
//          else if (string(i:i) .eq. '1') then
//             digit = 1
//          else if (string(i:i) .eq. '2') then
//             digit = 2
//          else if (string(i:i) .eq. '3') then
//             digit = 3
//          else if (string(i:i) .eq. '4') then
//             digit = 4
//          else if (string(i:i) .eq. '5') then
//             digit = 5
//          else if (string(i:i) .eq. '6') then
//             digit = 6
//          else if (string(i:i) .eq. '7') then
//             digit = 7
//          else if (string(i:i) .eq. '8') then
//             digit = 8
//          else if (string(i:i) .eq. '9') then
//             digit = 9
//          end if
//          number = number + digit * place(j)
//       end do
//       if (negative)  number = -number
//       return
//       end