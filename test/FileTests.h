///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileTests.h
 * \author tsharpe
 * \date Dec 6, 2011
 *
 * \brief
 */
#ifndef TEST_FILETESTS_H_
#define TEST_FILETESTS_H_

#include "system/file/Directory.h"
#include "system/file/File.h"
#include "system/file/SymLink.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>

#include <cxxtest/TestSuite.h>

class FileTests : public CxxTest::TestSuite
{
public:
    void testNameConstruction()
    {
        TS_ASSERT_EQUALS(File().toString(),"");
        TS_ASSERT_EQUALS(File("/").toString(),"/");
        TS_ASSERT_EQUALS(File("//").toString(),"/");
        TS_ASSERT_EQUALS(File("//./").toString(),"/");
        TS_ASSERT_EQUALS(File("//.").toString(),"/");

        File fooBar("foo/bar");
        TS_ASSERT_EQUALS(fooBar.toString(),"foo/bar");
        TS_ASSERT_EQUALS(File(std::string("foo/bar")).toString(),"foo/bar");

        File baz("baz");
        baz = fooBar;
        TS_ASSERT_EQUALS(baz.toString(),"foo/bar");

        File bah(fooBar);
        TS_ASSERT_EQUALS(bah.toString(),"foo/bar");

        TS_ASSERT_EQUALS(Directory().toString(),"/");
        TS_ASSERT_EQUALS(Directory("/").toString(),"/");
        TS_ASSERT_EQUALS(Directory(".").toString(),".");
        TS_ASSERT_EQUALS(Directory("foo/bar").toString(),"foo/bar");
    }

    void testStat()
    {
        Directory tmp("./tmp");
        File foo(tmp.file("foo"));
        TS_ASSERT_EQUALS(foo,"./tmp/foo");
        File bar(tmp.file("bar"));
        TS_ASSERT_EQUALS(bar,"./tmp/bar");
        SymLink baz(tmp.file("baz"));
        TS_ASSERT_EQUALS(baz,"./tmp/baz");

        tmp.create();
        TS_ASSERT(tmp.isValid());
        close(creat(foo.toString().c_str(),0777));
        time_t now = time(0);
        baz.setTarget("xyzzy");
        TS_ASSERT(!baz.isValid());
        TS_ASSERT(baz.isLink());

        File bareFoo("foo");
        baz.setTarget(bareFoo,true);
        TS_ASSERT(baz.isValid());

        TS_ASSERT_EQUALS(baz.target(),bareFoo);
        TS_ASSERT(baz.isSameFile(foo));

        TS_ASSERT(tmp.exists());
        TS_ASSERT(foo.exists());
        TS_ASSERT(!bar.exists());
        TS_ASSERT(baz.exists());

        TS_ASSERT(tmp.stat() != 0);
        TS_ASSERT(foo.stat() != 0);
        TS_ASSERT(bar.stat() == 0);

        TS_ASSERT_EQUALS(tmp.type(),File::DIRECTORY);
        TS_ASSERT_EQUALS(foo.type(),File::REGULAR);
        TS_ASSERT_EQUALS(bar.type(),File::NOT_A_FILE);
        TS_ASSERT_EQUALS(baz.type(),File::SYM_LINK);

        TS_ASSERT(tmp.isDir());
        TS_ASSERT(!foo.isDir());
        TS_ASSERT(!bar.isDir());
        TS_ASSERT(!baz.isDir());

        TS_ASSERT(!tmp.isLink());
        TS_ASSERT(!foo.isLink());
        TS_ASSERT(!bar.isLink());
        TS_ASSERT(baz.isLink());

        TS_ASSERT_DIFFERS(tmp.filesize(),0L);
        TS_ASSERT_EQUALS(foo.filesize(),0L);
        TS_ASSERT_EQUALS(bar.filesize(),-1L);
        TS_ASSERT_EQUALS(baz.filesize(),0L);

        TS_ASSERT_LESS_THAN(now-foo.modifyTime(),1L);

        foo.remove();
        TS_ASSERT_EQUALS(foo.type(),File::NOT_A_FILE);
        baz.remove();
        TS_ASSERT_EQUALS(baz.type(),File::NOT_A_FILE);
        tmp.remove();
        TS_ASSERT_EQUALS(tmp.type(),File::NOT_A_FILE);
    }

    void testPathManipulation()
    {
        Directory cwd(".");
        TS_ASSERT_EQUALS(cwd,".");
        Directory root("/");
        TS_ASSERT_EQUALS(root,"/");
        Directory dir("/home/user/boat/.");
        TS_ASSERT_EQUALS(dir,"/home/user/boat");
        File file(dir.file("file.txt"));
        TS_ASSERT_EQUALS(file,"/home/user/boat/file.txt");
        File hidden(".hidden");

        TS_ASSERT_EQUALS(cwd.filename(),".");
        TS_ASSERT_EQUALS(root.filename(),"");
        TS_ASSERT_EQUALS(dir.filename(),"boat");
        TS_ASSERT_EQUALS(file.filename(),"file.txt");
        TS_ASSERT_EQUALS(hidden.filename(),".hidden");

        TS_ASSERT_EQUALS(cwd.dirname(),".");
        TS_ASSERT_EQUALS(root.dirname(),"");
        TS_ASSERT_EQUALS(dir.dirname(),"/home/user");
        TS_ASSERT_EQUALS(file.dirname(),"/home/user/boat");
        TS_ASSERT_EQUALS(hidden.dirname(),".");

        TS_ASSERT_EQUALS(dir.extension(),"");
        TS_ASSERT_EQUALS(file.extension(),".txt");
        TS_ASSERT_EQUALS(hidden.extension(),"");

        TS_ASSERT_EQUALS(dir.removeExtension(),dir);
        TS_ASSERT_EQUALS(file.removeExtension(),"/home/user/boat/file");
        TS_ASSERT_EQUALS(hidden.removeExtension(),hidden);

        TS_ASSERT_EQUALS(dir.addExtension(".new"),"/home/user/boat.new");
        TS_ASSERT_EQUALS(file.addExtension(".new"),"/home/user/boat/file.txt.new");
        TS_ASSERT_EQUALS(hidden.addExtension(".new"),".hidden.new");

        TS_ASSERT_EQUALS(dir.changeExtension(".new"),"/home/user/boat.new");
        TS_ASSERT_EQUALS(file.changeExtension(".new"),"/home/user/boat/file.new");
        TS_ASSERT_EQUALS(hidden.changeExtension(".new"),".hidden.new");
    }
};

#endif /* TEST_FILETESTS_H_ */
